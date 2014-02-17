/*******************************************************************************
    Copyright (C) 2014 Brain Research Institute, Melbourne, Australia
    
    Permission is hereby granted under the Patent Licence Agreement between
    the BRI and Siemens AG from July 3rd, 2012, to Siemens AG obtaining a
    copy of this software and associated documentation files (the
    "Software"), to deal in the Software without restriction, including
    without limitation the rights to possess, use, develop, manufacture,
    import, offer for sale, market, sell, lease or otherwise distribute
    Products, and to permit persons to whom the Software is furnished to do
    so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included
    in all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
    OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*******************************************************************************/


#ifndef __mrtrix_thread_queue_h__
#define __mrtrix_thread_queue_h__

#include <stack>

#include "ptr.h"
#include "thread/condition.h"
#include "thread/exec.h"

#define MRTRIX_QUEUE_DEFAULT_CAPACITY 128
#define MRTRIX_QUEUE_DEFAULT_BATCH_SIZE 128

namespace MR
{
  namespace Thread
  {


    //* \cond skip
    namespace {

      /********************************************************************
       * convenience Functor classes for use in Thread::run_queue()
       ********************************************************************/

      template <class Functor> 
        class __Multi {
          public:
            __Multi (Functor& object, size_t number) : obj (object), num (number) { }
            template <class X> bool operator() (const X& x) { assert (0); return false; }
            template <class X> bool operator() (X& x) { assert (0); return false; }
            template <class X, class Y> bool operator() (const X& x, Y& y) { assert (0); return false; }
            Functor& obj;
            size_t num;
        };

      template <class Item> 
        class __Batch {
          public:
            __Batch (size_t number) : num (number) { }
            size_t num;
        };



      // to handle batched / unbatched seamlessly:
      template <class X> class __item { public: typedef X type; }; 
      template <class X> class __item < __Batch<X> > { public: typedef X type; };

      // to get multi/single job/functor seamlessly:
      template <class X>
        class __job
        {
          public:
            typedef X type;
            static X& functor (const X& job) { return const_cast<X&> (job); }
        };

      template <class X>
        class __job < __Multi<X> >
        {
          public:
            typedef X type;
            static X& functor (const __Multi<X>& job) { return const_cast<X&> (job.obj); }
        };

      // to launch multi/single job/functor seamlessly:
      template <class F, class T> 
        class __exec 
        { 
          public:
            __exec (const F& f, T& t, const char* name) : launcher (t, name) { }
          private:
            Exec launcher;
        };

      template <class F, class T> 
        class __exec <__Multi<F>, T>
        { 
          public:
            __exec (const __Multi<F>& f, T& t, const char* name) : 
              array (t, f.num),
              launcher (array, name) { }
          private:
            Array<T> array;
            Exec launcher;
        };

    }

    //! \endcond 




    /** \addtogroup Thread
     * @{ */



    /** \defgroup thread_queue Thread-safe queue
     * \brief Functionality for thread-safe parallel processing of queued items 
     *
     * These functions and classes provide functionality for one or more \e
     * source threads to feed items into a first-in first-out queue, and one or
     * more \e sink threads to consume items. This pipeline can also extend
     * over two queues, with one or more \e pipe threads consuming items of one
     * type from the first queue, and feeding items of another type onto the
     * second queue. 
     *
     * As a graphical representation of the pipeline, the following workflows
     * can be achieved:
     *
     * \code
     *     [source] \               / [sink]
     *     [source] -- queue<item> -- [sink]
     *     [source] /               \ [sink]
     *        ..                        ..
     *     N_source                   N_sink
     * \endcode
     *
     * or for a deeper pipeline:
     *
     * \code
     *     [source] \                / [pipe]  \                 / [sink]
     *     [source] -- queue<item1> -- [pipe]  -- queue<item2>  -- [sink]
     *     [source] /                \ [pipe]  /                 \ [sink]
     *        ..                         ..                          ..    
     *     N_source                    N_pipe                      N_sink
     * \endcode
     *
     * By default, items are push to and pulled from the queue one by one. In
     * situations where the amount of processing per item is small, items can
     * be sent in batches to reduce the overhead of thread management (mutex
     * locking/unlocking, etc). 
     *
     * The simplest way to use this functionality is via the
     * Thread::run_queue() and associated Thread::multi() and Thread::batch()
     * functions. In complex situations, it may be necessary to use the
     * Thread::Queue class directly, although that should very rarely (if ever)
     * be needed. 
     *
     * \sa Thread::run_queue()
     * \sa Thread::Queue
     *
     * @{ */



    //! A first-in first-out thread-safe item queue
    /*! This class implements a thread-safe means of pushing data items into a
     * queue, so that they can each be processed in one or more separate
     * threads. 
     *
     * \note In practice, it is almost always simpler to use the convenience
     * function Thread::run_queue(). You should never need to use the
     * Thread::Queue directly unless you have a very unusual situation.
     *
     * \section thread_queue_usage Usage overview
     *
     * Thread::Queue has somewhat unusual usage, which consists of the following
     * steps:
     * - Create an instance of a Thread::Queue
     * - Create one or more instances of the corresponding
     * Thread::Queue::Writer class, each constructed with a reference to the
     * queue. Each of these instances will automatically notify the queue that
     * its corresponding thread will be writing to the queue.
     * - Create one or more instances of the corresponding
     * Thread::Queue::Reader class, each constructed with a reference to the
     * queue. Each of these instances will automatically notify the queue that
     * its corresponding thread will be reading from the queue.
     * - Launch all threads, one per instance of Thread::Queue::Writer or
     *   Thread::Queue::Reader. Note that one of these threads can be the
     *   current thread - simply invoke the respective functor's execute()
     *   method directly once all other threads have been launched.
     * - Within the execute() method of each thread with a
     * Thread::Queue::Writer:
     *   - create an instance of Thread::Queue::Writer::Item, constructed from
     *   the corresponding Thread::Queue::Writer;
     *   - perform processing in a loop:
     *     - prepare the item using pointer semantics (i.e. *item or
     *     item->method());
     *     - use the write() method of this class to write to the queue;
     *     - break out of loop if write() returns \c false.
     *   - when the execute() method returns, the destructor of the
     *   Thread::Queue::Writer::Item class will notify the queue that its
     *   thread has finished writing to the queue.
     * - Within the execute() method of each thread with a
     * Thread::Queue::Reader:
     *   - create an instance of Thread::Queue::Reader::Item, constructed from
     *   the corresponding Thread::Queue::Reader;
     *   - perform processing in a loop:
     *     - use the read() method of this class to read the next item from the
     *     queue;
     *     - break out of the loop if read() returns \c false;
     *     - process the item using pointer semantics (i.e. *item or
     *     item->method()).
     *   - when the execute() method returns, the destructor of the
     *   Thread::Queue::Reader::Item class will notify the queue that its
     *   thread has finished reading from the queue.
     * - If all reader threads have returned, the queue will notify all writer
     * threads that processing should stop, by returning \c false from the next
     * write attempt.
     * - If all writer threads have returned and no items remain in the queue,
     * the queue will notify all reader threads that processing should stop, by
     * returning \c false from the next read attempt.
     *
     * The additional member classes are designed to be used in conjunction
     * with the MRtrix multi-threading interface. In this system, each thread
     * corresponds to an instance of a functor class, and its execute() method
     * is the function that will be run within the thread (see Thread::Exec for
     * details). For this reason:
     * - The Thread::Queue instance is designed to be created before any of the
     *   threads.
     * - The Thread::Queue::Writer and Thread::Queue::Reader classes are
     *   designed to be used as members of each functor, so that each functor
     *   must construct these classes from a reference to the queue within
     *   their own constructor. This ensures each thread registers their
     *   intention to read or write with the queue \e before their thread is
     *   launched.
     * - The Thread::Queue::Writer::Item and Thread::Queue::Reader::Item
     *   classes are designed to be instantiated within each functor's
     *   execute() method. They must be constructed from a reference to a
     *   Thread::Queue::Writer or Thread::Queue::Reader respectively, ensuring
     *   no reads or write can take place without having registered with the
     *   queue. Their destructors will also unregister from the queue, ensuring
     *   that each thread unregisters as soon as the execute() method returns,
     *   and hence \e before the thread exits.
     *
     * The Queue class performs all memory management for the items in the
     * queue. For this reason, the items are accessed via the Writer::Item &
     * Reader::Item classes. This allows items to be recycled once they have
     * been processed, reducing overheads associated with memory
     * allocation/deallocation. 
     *
     * \note It is important that all instances of Thread::Queue::Writer and
     * Thread::Queue::Reader are created \e before any of the threads are
     * launched, to avoid any race conditions at startup.
     *
     * The use of Thread::Queue is best illustrated with an example:
     * \code
     * // the type of objects that will be sent through the queue:
     * class Item {
     *   public:
     *     ...
     *     // data members
     *     ...
     * };
     *
     *
     * // The use a typedef is recommended to help with readability (and typing!):
     * typedef Thread::Queue<Item> MyQueue;
     *
     *
     * // this class will write to the queue:
     * class Sender {
     *   public:
     *     // construct the 'writer' member in the constructor:
     *     Sender (MyQueue& queue) : writer (queue) { }
     *
     *     void execute () {
     *       // use a local instance of Thread::Queue<Item>::Writer::Item to write to the queue:
     *       MyQueue::Writer::Item item (writer);
     *       while (need_more_items()) {
     *         ...
     *         // prepare item
     *         *item = something();
     *         item->set (something_else);
     *         ...
     *         if (!item.write()) break; // break if write() returns false
     *       }
     *     }
     *
     *   private:
     *     MyQueue::Writer writer;
     * };
     *
     *
     * // this class will read from the queue:
     * class Receiver {
     *   public:
     *     // construct the 'reader' member in the constructor:
     *     Receiver (MyQueue& queue) : reader (queue) { }
     *
     *     void execute () {
     *       // use a local instance of Thread::Queue<Item>::Reader::Item to read from the queue:
     *       MyQueue::Reader::Item item (reader);
     *       while ((item.read())) { // break when read() returns false
     *         ...
     *         // process item
     *         do_something (*item);
     *         if (item->status()) report_error();
     *         ...
     *         if (enough_items()) return;
     *       }
     *     }
     *
     *   private:
     *     MyQueue::Reader reader;
     * };
     *
     *
     * // this is where the queue and threads are created:
     * void my_function () {
     *   // create an instance of the queue:
     *   MyQueue queue;
     *
     *   // create all functors from a reference to the queue:
     *   Sender sender (queue);
     *   Receiver receiver (queue);
     *
     *   // once all functors are created, launch their corresponding threads:
     *   Thread::Exec sender_thread (sender);
     *   Thread::Exec receiver_thread (receiver);
     * }
     * \endcode
     *
     * \section thread_queue_rationale Rationale for the Writer, Reader, and Item member classes
     *
     * The motivation for the use of additional member classes to perform the
     * actual process of writing and reading to and from the queue is related
     * to the need to keep track of the number of processes currently using the
     * queue. This is essential to ensure that threads are notified when the
     * queue is closed. This happens either when all readers have finished
     * reading; or when all writers have finished writing and no items are left
     * in the queue. This is complicated by the need to ensure that the various
     * operations are called in the right order to avoid deadlocks.
     *
     * There are essentially 4 operations that need to take place:
     * - registering an intention to read/write from/to the queue
     * - launching the corresponding thread
     * - unregistering from the queue
     * - terminating the thread
     *
     * For proper multi-threaded operations, these operations must take place
     * in the order above. Moreover, each operation must be completed for
     * all users of the queue before any of them can perform the next
     * operation. The use of additional member classes ensures that threads
     * have to register their intention to read or write from the queue, and
     * that they unregister from the queue once their processing is done.
     *
     * While this could have been achieved simply with the appropriate member
     * functions (i.e. register(), unregister(), %read() & %write() methods in
     * the main Queue class), this places a huge burden on the developer to get
     * it right. Using these member functions reduces the chance of coding
     * errors, and in fact reduces the total amount of code that needs to be
     * written to use the Queue in a safe manner.
     *
     * The Item classes additionally simplify the memory management of the
     * items in the queue, by preventing direct access to the underlying
     * pointers, and ensuring the Queue itself is responsible for all
     * allocation and deallocation of items as needed.
     *
     * \sa Thread::run_queue()
     */
    template <class T> class Queue
    {
      public:
        //! Construct a Queue of items of type \c T
        /*! \param description a string identifying the queue for degugging purposes
         * \param buffer_size the maximum number of items that can be pushed onto the queue before
         * blocking. If a thread attempts to push more data onto the queue when the
         * queue already contains this number of items, the thread will block until
         * at least one item has been popped.  By default, the buffer size is
         * MRTRIX_QUEUE_DEFAULT_CAPACITY items.
         */
        Queue (const std::string& description = "unnamed", size_t buffer_size = MRTRIX_QUEUE_DEFAULT_CAPACITY) :
          more_data (mutex),
          more_space (mutex),
          buffer (new T* [buffer_size]),
          front (buffer),
          back (buffer),
          capacity (buffer_size),
          writer_count (0),
          reader_count (0),
          name (description) {
          assert (capacity > 0);
        }

        //! needed for Thread::run_queue()
        Queue (const T& item_type, const std::string& description = "unnamed", size_t buffer_size = MRTRIX_QUEUE_DEFAULT_CAPACITY) :
          more_data (mutex),
          more_space (mutex),
          buffer (new T* [buffer_size]),
          front (buffer),
          back (buffer),
          capacity (buffer_size),
          writer_count (0),
          reader_count (0),
          name (description) {
          assert (capacity > 0);
        }


        ~Queue () {
          delete [] buffer;
        }

        //! This class is used to register a writer with the queue
        /*! Items cannot be written directly onto a Thread::Queue queue. An
         * object of this class must first be instanciated to notify the queue
         * that a section of code will be writing to the queue. The actual
         * process of writing items to the queue is done via the Writer::Item
         * class.
         *
         * \sa Thread::Queue for more detailed information and examples.
         * \sa Thread::run_queue() for a much more user-friendly way of setting
         * up a queue.  */
        class Writer
        {
          public:
            //! Register a Writer object with the queue
            /*! The Writer object will register itself with the queue as a
             * writer. */
            Writer (Queue<T>& queue) : Q (queue) {
              Q.register_writer();
            }
            Writer (const Writer& W) : Q (W.Q) {
              Q.register_writer();
            }

            //! This class is used to write items to the queue
            /*! Items cannot be written directly onto a Thread::Queue queue. An
             * object of this class must be instantiated and used to write to the
             * queue.
             *
             * \sa Thread::Queue for more detailed information and examples.
             * \sa Thread::run_queue() for a much more user-friendly way of setting
             * up a queue.  */
            class Item
            {
              public:
                //! Construct a Writer::Item object
                /*! The Writer::Item object can only be instantiated from a
                 * Writer object, ensuring that the corresponding section of code
                 * has already registered as a writer with the queue. The
                 * destructor for this object will unregister from the queue.
                 *
                 * \note There should only be one Writer::Item object per Writer.
                 * */
                Item (const Writer& writer) : Q (writer.Q), p (Q.get_item()) { }
                //! Unregister the parent Writer from the queue
                ~Item () {
                  Q.unregister_writer();
                }
                //! Push the item onto the queue
                bool write () {
                  return (Q.push (p));
                }
                T& operator*() const throw ()   {
                  return (*p);
                }
                T* operator->() const throw ()  {
                  return (p);
                }
              private:
                Queue<T>& Q;
                T* p;
            };

          private:
            Queue<T>& Q;
        };


        //! This class is used to register a reader with the queue
        /*! Items cannot be read directly from a Thread::Queue queue. An
         * object of this class must be instanciated to notify the queue
         * that a section of code will be reading from the queue. The actual
         * process of reading items from the queue is done via the Reader::Item
         * class.
         *
         * \sa Thread::Queue for more detailed information and examples.
         * \sa Thread::run_queue() for a much more user-friendly way of setting
         * up a queue.  */
        class Reader
        {
          public:
            //! Register a Reader object with the queue.
            /*! The Reader object will register itself with the queue as a
             * reader. */
            Reader (Queue<T>& queue) : Q (queue) {
              Q.register_reader();
            }
            Reader (const Reader& reader) : Q (reader.Q) {
              Q.register_reader();
            }

            //! This class is used to read items from the queue
            /*! Items cannot be read directly from a Thread::Queue queue. An
             * object of this class must be instanciated and used to read from the
             * queue.
             *
             * \sa Thread::Queue for more detailed information and examples.
             * \sa Thread::run_queue() for a much more user-friendly way of setting
             * up a queue.  */
            class Item
            {
              public:
                //! Construct a Reader::Item object
                /*! The Reader::Item object can only be instantiated from a
                 * Reader object, ensuring that the corresponding section of code
                 * has already registered as a reader with the queue. The
                 * destructor for this object will unregister from the queue.
                 *
                 * \note There should only be one Reader::Item object per
                 * Reader. */
                Item (const Reader& reader) : Q (reader.Q), p (NULL) { }
                //! Unregister the parent Reader from the queue
                ~Item () {
                  Q.unregister_reader();
                }
                //! Get next item from the queue
                bool read () {
                  return (Q.pop (p));
                }
                T& operator*() const throw ()   {
                  return (*p);
                }
                T* operator->() const throw ()  {
                  return (p);
                }
                bool operator! () const throw () {
                  return !p;
                }
              private:
                Queue<T>& Q;
                T* p;
            };
          private:
            Queue<T>& Q;
        };

        //! Print out a status report for debugging purposes
        void status () {
          Mutex::Lock lock (mutex);
          std::cerr << "Thread::Queue \"" + name + "\": "
                    << writer_count << " writer" << (writer_count > 1 ? "s" : "") << ", "
                    << reader_count << " reader" << (reader_count > 1 ? "s" : "") << ", items waiting: " << size() << "\n";
        }


      private:
        Mutex mutex;
        Cond more_data, more_space;
        T** buffer;
        T** front;
        T** back;
        size_t capacity;
        size_t writer_count, reader_count;
        std::stack<T*,std::vector<T*> > item_stack;
        VecPtr<T> items;
        std::string name;

        Queue (const Queue& queue) {
          assert (0);
        }
        Queue& operator= (const Queue& queue) {
          assert (0);
          return (*this);
        }

        void register_writer ()   {
          Mutex::Lock lock (mutex);
          ++writer_count;
        }
        void unregister_writer () {
          Mutex::Lock lock (mutex);
          assert (writer_count);
          --writer_count;
          if (!writer_count) {
            DEBUG ("no writers left on queue \"" + name + "\"");
            more_data.broadcast();
          }
        }
        void register_reader ()   {
          Mutex::Lock lock (mutex);
          ++reader_count;
        }
        void unregister_reader () {
          Mutex::Lock lock (mutex);
          assert (reader_count);
          --reader_count;
          if (!reader_count) {
            DEBUG ("no readers left on queue \"" + name + "\"");
            more_space.broadcast();
          }
        }

        bool empty () const {
          return (front == back);
        }
        bool full () const {
          return (inc (back) == front);
        }
        size_t size () const {
          return ( (back < front ? back+capacity : back) - front);
        }

        T* get_item () {
          Mutex::Lock lock (mutex);
          T* item = new T;
          items.push_back (item);
          return (item);
        }

        bool push (T*& item) {
          Mutex::Lock lock (mutex);
          while (full() && reader_count) more_space.wait();
          if (!reader_count) return (false);
          *back = item;
          back = inc (back);
          more_data.signal();
          if (item_stack.empty()) {
            item = new T;
            items.push_back (item);
          }
          else {
            item = item_stack.top();
            item_stack.pop();
          }
          return (true);
        }

        bool pop (T*& item) {
          Mutex::Lock lock (mutex);
          if (item) 
            item_stack.push (item);
          item = NULL;
          while (empty() && writer_count) more_data.wait();
          if (empty() && !writer_count) return (false);
          item = *front;
          front = inc (front);
          more_space.signal();
          return (true);
        }

        T** inc (T** p) const {
          ++p;
          if (p >= buffer + capacity) p = buffer;
          return (p);
        }
    };






     //* \cond skip

    template <class T> class Queue< __Batch<T> >
    {
      private:
        typedef std::vector<T> BatchType;
        typedef Queue<BatchType> BatchQueue;

      public:
        Queue (const __Batch<T>& item_type, const std::string& description = "unnamed", size_t buffer_size = MRTRIX_QUEUE_DEFAULT_CAPACITY) :
          batch_queue (description, buffer_size),
          batch_size (item_type.num) { }


        class Writer
        {
          public:
            Writer (Queue< __Batch<T> >& queue) : batch_writer (queue.batch_queue), batch_size (queue.batch_size) { }

            class Item
            {
              public:
                Item (const Writer& writer) : 
                  batch_item (writer.batch_writer), batch_size (writer.batch_size), n (0) { 
                    batch_item->resize (batch_size);
                }
                ~Item () {
                  if (n) {
                    batch_item->resize (n);
                    batch_item.write();
                  }
                }
                bool write () {
                  if (++n >= batch_size) {
                    if (!batch_item.write()) 
                      return false;
                    n = 0;
                    batch_item->resize (batch_size);
                  }
                  return true;
                }
                T& operator*() const throw ()   {
                  return (*batch_item)[n];
                }
                T* operator->() const throw ()  {
                  return &((*batch_item)[n]);
                }
              private:
                typename BatchQueue::Writer::Item batch_item;
                const size_t batch_size;
                size_t n;
            };

          private:
            typename BatchQueue::Writer batch_writer;
            const size_t batch_size;
        };


        class Reader
        {
          public:
            Reader (Queue< __Batch<T> >& queue) : batch_reader (queue.batch_queue), batch_size (queue.batch_size) { }

            class Item
            {
              public:
                Item (const Reader& reader) : batch_item (reader.batch_reader), batch_size (reader.batch_size), n (0) { }
                bool read () {
                  if (!batch_item) 
                    return batch_item.read();

                  if (++n >= batch_item->size()) {
                    if (!batch_item.read()) 
                      return false;
                    n = 0;
                  }
                  return true;
                }
                T& operator*() const throw ()   {
                  return (*batch_item)[n];
                }
                T* operator->() const throw ()  {
                  return &((*batch_item)[n]);
                }
              private:
                typename BatchQueue::Reader::Item batch_item;
                const size_t batch_size;
                size_t n;
            };
          private:
            typename BatchQueue::Reader batch_reader;
            const size_t batch_size;
        };

        void status () { batch_queue.status(); }


      private:
        BatchQueue batch_queue;
        const size_t batch_size;
    };







    /*! wrapper classes to extend simple functors designed for use with
     * Thread::run_queue with functionality needed for use with Thread::Queue */
    namespace {

       template <class Type, class Functor>
         class __Source
         {
           public:
             __Source (Queue<Type>& queue, Functor& functor) :
               writer (queue), func (__job<Functor>::functor (functor)) { }
             __Source (const __Source& S) :
               writer (S.writer), funcp (new Functor (S.func)), func (*funcp) { }

             void execute () {
               typename Queue<Type>::Writer::Item out (writer);
               do {
                 if (!func (*out)) 
                   return;
               } while (out.write());
             }

           private:
             typename Queue<Type>::Writer writer;
             Ptr<typename __job<Functor>::type> funcp;
             typename __job<Functor>::type& func;
         };


       template <class Type1, class Functor, class Type2>
         class __Pipe
         {
           public:
             typedef typename __job<Functor>::type Functype;

             __Pipe (Queue<Type1>& queue_in, Functor& functor, Queue<Type2>& queue_out) :
               reader (queue_in), writer (queue_out), func (__job<Functor>::functor (functor)) { }
             __Pipe (const __Pipe& P) :
               reader (P.reader), writer (P.writer), funcp (new Functor (P.func)), func (*funcp) { }

             void execute () {
               typename Queue<Type1>::Reader::Item in (reader);
               typename Queue<Type2>::Writer::Item out (writer);
               do {
                 if (!in.read()) return;
                 if (!func (*in, *out))
                   return;
               }
               while (out.write());
             }

           private:
             typename Queue<Type1>::Reader reader;
             typename Queue<Type2>::Writer writer;
             Ptr<typename __job<Functor>::type> funcp;
             typename __job<Functor>::type& func;
         };



       template <class Type, class Functor>
         class __Sink
         {
           public:
             typedef typename __job<Functor>::type Functype;

             __Sink (Queue<Type>& queue, Functor& functor) :
               reader (queue), func (__job<Functor>::functor (functor)) { }
             __Sink (const __Sink& S) :
               reader (S.reader), funcp (new Functor (S.func)), func (*funcp) { }

             void execute () {
               typename Queue<Type>::Reader::Item in (reader);
               while (in.read()) {
                 if (!func (*in))
                   return;
               }
             }

           private:
             typename Queue<Type>::Reader reader;
             Ptr<typename __job<Functor>::type> funcp;
             typename __job<Functor>::type& func;
         };


    }


    //! \endcond








    //! used to request multiple threads of the corresponding functor
    /*! This function is used in combination with Thread::run_queue to request
     * that the functor \a object be run in parallel using \a number threads of
     * execution (defaults to Thread::number_of_threads()). 
     * \sa Thread::run_queue() */
    template <class Functor>
      inline __Multi<Functor> multi (Functor& object, size_t number = number_of_threads()) 
      {
        return __Multi<Functor> (object, number);
      }

 
    //! used to request batched processing of items
    /*! This function is used in combination with Thread::run_queue to request
     * that the items \a object be processed in batches of \a number items 
     * (defaults to MRTRIX_QUEUE_DEFAULT_BATCH_SIZE).
     * \sa Thread::run_queue() */
    template <class Item>
      inline __Batch<Item> batch (const Item& object, size_t number = MRTRIX_QUEUE_DEFAULT_BATCH_SIZE) 
      {
        return __Batch<Item> (number);
      }






    //! convenience function to set up and run a 2-stage multi-threaded pipeline.
    /*! This function (and its 3-stage equivalent 
     * Thread::run_queue(const Source&, const Type1&, const Pipe&, const Type2&, const Sink&, size_t)) 
     * simplify the process of setting up a multi-threaded processing chain
     * that should meet most users' needs.
     *
     * The arguments to this function correspond to an instance of the Source, 
     * the Sink, and optionally the Pipe functors, in addition to an instance
     * of the Items to be passed through each stage of the pipeline - these are
     * provided purely to specify the type of object to pass through the
     * queue(s).
     *
     * \note While the functors are passed by const-reference, you  should not
     * rely on them being unmodified - they will immediately be cast back
     * to their non-const version within Thread::run_queue(). The reason for
     * the const declaration is to support passing the lightweight temporary
     * structures produced by Thread::multi(), which have to be passed by
     * const-reference as mandated by the C++ standard.
     *
     * \section thread_run_queue_functors Functors
     *
     * The 3 types of functors each have a specific purpose, and corresponding
     * requirements as described below:
     *
     * \par Source: the input functor 
     * The Source class must at least provide the method:
     * \code 
     * bool operator() (Type& item);
     * \endcode
     * This function prepares the \a item passed to it, and should return \c
     * true if further items need to be processed, or \c false to signal that
     * no further items are to be sent through the queue (at which point the
     * corresponding thread(s) will exit).
     *
     * \par Sink: the output functor 
     * The Sink class must at least provide the method:
     * \code 
     * bool operator() (const Type& item);
     * \endcode
     * This function processes the \a item passed to it, and should return \c
     * true when ready to process further items, or \c false to signal the end
     * of processing (at which point the corresponding thread(s) will exit).
     *
     * \par Pipe: the processing functor (for 3-stage pipeline only)
     * The Pipe class must at least provide the method:
     * \code 
     * bool operator() (const Type1& item_in, Type2& item_out);
     * \endcode
     * This function processes the \a item_in passed to it, and prepares
     * \a item_out for the next stage of the pipeline. It should return \c
     * true when ready to process further items, or \c false to signal the end
     * of processing (at which point the corresponding thread(s) will exit).
     *
     * \section thread_run_queue_example Simple example
     *
     * This is a simple demo application that generates a linear sequence of
     * numbers and sums them up:
     *
     * \code
     * const size_t max_count;
     *
     * // the functor that will generate the items:
     * class Source {
     *   public:
     *     Source () : count (0) { }
     *     bool operator() (size_t& item) {
     *       item.value = count++;
     *       return count < max_count; // stop when max_count is reached
     *     }
     * };
     *
     * // the functor that will consume the items:
     * class Sink {
     *   public:
     *     Sink (size_t& total) : 
     *       grand_total (grand_total),
     *       total (0) { }
     *     ~Sink () { // update grand_total in destructor
     *       grand_total += total;
     *     }
     *     bool operator() (const size_t& item) {
     *       total += item.value();
     *       return true;
     *    }
     *  protected:
     *    size_t& grand_total;
     * };
     *
     * void run () 
     * {
     *   size_t grand_total = 0;
     *   Source source;
     *   Sink sink (grand_total);
     *
     *   // run a single-source => single-sink pipeline:
     *   Thread::run_queue (source, size_t(), sink);
     * }
     * \endcode
     *
     * \section thread_run_queue_multi Parallel execution of functors
     *
     * If a functor is to be run over multiple parallel threads of execution,
     * it should be wrapped in a call to Thread::multi() before being passed
     * to the Thread::run_queue() functions.  The Thread::run_queue() functions
     * will then create additional instances of the relevant functor using its
     * copy constructor, via the Thread::Array class; care should therefore be
     * taken to ensure that the functor's copy constructor behaves
     * appropriately.
     *
     * For example, using the code above:
     *
     * \code 
     * ...
     *
     * void run () 
     * {
     *   ...
     *
     *   // run a single-source => multi-sink pipeline:
     *   Thread::run_queue (source, size_t(), Thread::multi (sink));
     * }
     * \endcode
     *
     * For the functor that is being multi-threaded, the default number of
     * threads instantiated will depend on the "NumberOfThreads" entry in the
     * MRtrix confugration file, or can be set at the command-line using the
     * -nthreads option. This number can also be set as additional optional
     * argument to Thread::multi().
     *
     * Note that any functor can be parallelised in this way. In the example
     * above, the Source functor could have been wrapped in Thread::multi()
     * instead if this was the behaviour required:
     *
     * \code 
     * ...
     *
     * void run () 
     * {
     *   ...
     *
     *   // run a multi-source => single-sink pipeline:
     *   Thread::run_queue (Thread::multi (source), size_t(), sink);
     * }
     * \endcode
     *
     *
     * \section thread_run_queue_batch Batching items
     *
     * In cases where the amount of processing per item is small, the overhead
     * of managing the concurrent access to the various queues from all the
     * threads may become prohibitive (see \ref multithreading for details). In
     * this case, it is a good idea to process the items in batches, which
     * drastically reduces the number of accesses to the queue. This can be
     * done by wrapping the items in a call to Thread::batch():
     *
     * \code 
     * ...
     *
     * void run () 
     * {
     *   ...
     *
     *   // run a single-source => multi-sink pipeline on batches of size_t items:
     *   Thread::run_queue (source, Thread::batch (size_t()), Thread::multi (sink));
     * }
     * \endcode
     *
     * By default, batches consist of MRTRIX_QUEUE_DEFAULT_BATCH_SIZE items
     * (defined as 128). This can be set explicitly by providing the desired
     * size as an additional argument to Thread::batch():
     *
     * \code 
     * ...
     *
     * void run () 
     * {
     *   ...
     *
     *   // run a single-source => multi-sink pipeline on batches of 1024 size_t items:
     *   Thread::run_queue (source, Thread::batch (size_t(), 1024), Thread::multi (sink));
     * }
     * \endcode
     *
     * Obviously, Thread::multi() and Thread::batch() can be used in any
     * combination to perform the operations required. 
     */

    template <class Source, class Type, class Sink>
      inline void run_queue (
          const Source& source, 
          const Type& item_type, 
          const Sink& sink, 
          size_t capacity = MRTRIX_QUEUE_DEFAULT_CAPACITY)
      {
        if (number_of_threads() == 0) {
          typename __item<Type>::type item;
          while (const_cast<Source&> (source) (item)) 
            if (!const_cast<Sink&> (sink) (item))
              return;
          return;
        }

        typedef typename __job<Source>::type SourceType;
        typedef typename __job<Sink>::type SinkType;

         Queue<Type> queue (item_type, "source->sink", capacity);
         __Source<Type,SourceType> source_functor (queue, __job<Source>::functor (source));
         __Sink<Type,SinkType>     sink_functor   (queue, __job<Sink>::functor (sink));

         __exec<Source, __Source<Type,SourceType> > source_threads (source, source_functor, "source");
         __exec<Sink, __Sink<Type,SinkType> >       sink_threads   (sink, sink_functor, "sink");
      }



    //! convenience functions to set up and run a 3-stage multi-threaded pipeline.
    /*! This function extends the 2-stage Thread::run_queue() function to allow
     * a 3-stage pipeline. For example, using the example from 
     * Thread::run_queue(), the following would add an additional stage to the
     * pipeline to double the numbers as they come through:
     *
     * \code
     *
     * ...
     *
     * class Pipe {
     *   public:
     *     bool operator() (const Item& item_in, Item& item_out) {
     *       item_out.n = 2 * item_in.n;
     *       return true;
     *     }
     * };
     *
     * ...
     *
     * void run () 
     * {
     *   ...
     *
     *   // run a single-source => multi-pipe => single-sink pipeline on batches of size_t items:
     *   Thread::run_queue (
     *       source, 
     *       Thread::batch (size_t()), 
     *       Thread::multi (sink)
     *       Thread::batch (size_t()), 
     *       sink);
     * }
     * \endcode
     *
     * Note that any functor can be executed in parallel (i.e. wrapped in
     * Thread::multi()), Items do not need to be of the same type, and can be
     * batched independently with any desired size. 
     * */
    template <class Source, class Type1, class Pipe, class Type2, class Sink>
      inline void run_queue (
          const Source& source,
          const Type1& item_type1, 
          const Pipe& pipe, 
          const Type2& item_type2, 
          const Sink& sink, 
          size_t capacity = MRTRIX_QUEUE_DEFAULT_CAPACITY)
      {
        if (number_of_threads() == 0) {
          typename __item<Type1>::type item1;
          typename __item<Type2>::type item2;
          while (const_cast <Source&> (source) (item1)) {
            if (!const_cast<Pipe&> (pipe) (item1, item2))
              return;
            if (!const_cast<Sink&> (sink) (item2))
              return;
          }
          return;
        }


        typedef typename __job<Source>::type SourceType;
        typedef typename __job<Pipe>::type PipeType;
        typedef typename __job<Sink>::type SinkType;

         Queue<Type1> queue1 (item_type1, "source->pipe", capacity);
         Queue<Type2> queue2 (item_type2, "pipe->sink", capacity);

         __Source<Type1,SourceType>   source_functor (queue1, __job<Source>::functor (source));
         __Pipe<Type1,PipeType,Type2> pipe_functor   (queue1, __job<Pipe>::functor (pipe), queue2);
         __Sink<Type2,SinkType>       sink_functor   (queue2, __job<Sink>::functor (sink));

         __exec<Source, __Source<Type1,SourceType> >   source_threads (source, source_functor, "source");
         __exec<Pipe, __Pipe<Type1,PipeType,Type2> >   pipe_threads (pipe, pipe_functor, "pipe");
         __exec<Sink, __Sink<Type2,SinkType> >         sink_threads (sink, sink_functor, "sink");
      }


    /** @} */
    /** @} */
  }
}

#endif


