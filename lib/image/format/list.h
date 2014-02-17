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


#ifndef __image_format_list_h__
#define __image_format_list_h__

#include "image/header.h"

#define DECLARE_IMAGEFORMAT(format, desc) \
  class format : public Base { \
    public:  \
      format () : Base (desc) { } \
      virtual RefPtr<Handler::Base> read (Header& H) const; \
      virtual bool check (Header& H, size_t num_axes) const; \
      virtual RefPtr<Handler::Base> create (Header& H) const; \
  }

namespace MR
{
  namespace Image
  {

    /*! \todo add Image::Format::Folder class to handle both DICOM, GE I-dot
     * and potentially Siemens Numaris 3 images. */

    //! Classes responsible for handling of specific image formats
    namespace Format
    {

      //! the interface for classes that support the various image formats.
      /*! All image formats supported by %MRtrix are handled by a class derived
       * from the Image::Format::Base. An instance of each of these classes is
       * added to the list in the file list.cpp. */
      class Base
      {
        public:
          Base (const char* desc) : description (desc) { }
          virtual ~Base() { }

          const char*  description; //!< a short human-readable description of the image format

          /*! \brief read image file(s) and fill the Image::Header \c H with the
           * appropriate information.
           *
           * This function will check whether this handler can read images in
           * the format suggested by the filename. It will then attempt to read
           * the corresponding image, and update the Image::Header \c H with
           * the relevant information.
           *
           * \return true if this Format handler can read this type of file,
           * false otherwise.  If true, this function should fill the Header \c
           * H with all the relevant information as read from the images before
           * returning.
           * \note this function should throw an Exception in case of error. */
          virtual RefPtr<Handler::Base> read (Header& H) const = 0;

          /*! \brief check whether the Image::Header \c H can be created using
           * this handler.
           *
           * This function will check whether this handler can write images in
           * the format suggested by the filename. It will then check whether
           * the format can handle the number of dimensions requested, and
           * modify the header appropriately if needed.
           *
           * \return true if this Format handler can write this type of file,
           * false otherwise.
           * \note this function should throw an Exception in case of error, or
           * if this image format cannot support the header information. */
          virtual bool check (Header& H, size_t num_axes) const = 0;

          /*! \brief create the image corresponding to the Image::Header \c H.
           *
           * This function will create images in the corresponding format,
           * assuming the header has been validated using the check() function
           * beforehand.
           *
           * \note this function should throw an Exception in case of error. */
          virtual RefPtr<Handler::Base> create (Header& H) const = 0;
      };

      DECLARE_IMAGEFORMAT (Pipe, "Internal pipe");
      DECLARE_IMAGEFORMAT (DICOM, "DICOM");
      DECLARE_IMAGEFORMAT (MRtrix, "MRtrix");
      DECLARE_IMAGEFORMAT (MRtrix_GZ, "MRtrix (GZip compressed)");
      DECLARE_IMAGEFORMAT (NIfTI, "NIfTI-1.1");
      DECLARE_IMAGEFORMAT (NIfTI_GZ, "NIfTI-1.1 (GZip compressed)");
      DECLARE_IMAGEFORMAT (Analyse, "AnalyseAVW / NIfTI-1.1");
      DECLARE_IMAGEFORMAT (MRI, "MRTools (legacy format)");
      DECLARE_IMAGEFORMAT (XDS, "XDS");
      DECLARE_IMAGEFORMAT (MGH, "MGH");
      DECLARE_IMAGEFORMAT (MGZ, "MGZ (compressed MGH)");

      /*! a list of all extensions for image formats that %MRtrix can handle. */
      extern const char* known_extensions[];

      /*! a list of all handlers for supported image formats. */
      extern const Base* handlers[];



      /*! \fn bool Analyse::read (Header& H) const
        \todo need to update Analyse Format handler to new NIfTI utils
       * framework. */

      /*! \fn void Analyse::create (Header& H) const
        \todo need to implement writing in Analyse Format (using the new
       * NIfTI utils framework. */

      //! \cond skip
      extern Format::MRtrix mrtrix_handler;
      //! \endcond
    }
  }
}



#endif
