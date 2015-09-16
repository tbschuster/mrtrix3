/*
   Copyright 2009 Brain Research Institute, Melbourne, Australia

   Written by J-Donald Tournier, 2014.

   This file is part of MRtrix.

   MRtrix is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   MRtrix is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

*/

//#define GL_DEBUG

#include <string>

#include "point.h"
#include "image/buffer.h"
#include "math/versor.h"

#include "gui/cursor.h"
#include "gui/projection.h"
#include "gui/dialog/file.h"
#include "gui/mrview/tool/roi_editor/roi.h"


namespace MR
{
  namespace GUI
  {
    namespace MRView
    {
      namespace Tool
      {

            


        ROI::ROI (Dock* parent) :
            Base (parent),
            in_insert_mode (false)
        {

          VBoxLayout* main_box = new VBoxLayout (this);
          HBoxLayout* layout = new HBoxLayout;
          layout->setContentsMargins (0, 0, 0, 0);
          layout->setSpacing (0);

          QPushButton* button = new QPushButton (this);
          button->setToolTip (tr ("New ROI"));
          button->setIcon (QIcon (":/new.svg"));
          connect (button, SIGNAL (clicked()), this, SLOT (new_slot ()));
          layout->addWidget (button, 1);

          button = new QPushButton (this);
          button->setToolTip (tr ("Open ROI"));
          button->setIcon (QIcon (":/open.svg"));
          connect (button, SIGNAL (clicked()), this, SLOT (open_slot ()));
          layout->addWidget (button, 1);

          save_button = new QPushButton (this);
          save_button->setToolTip (tr ("Save ROI"));
          save_button->setIcon (QIcon (":/save.svg"));
          save_button->setEnabled (false);
          connect (save_button, SIGNAL (clicked()), this, SLOT (save_slot ()));
          layout->addWidget (save_button, 1);

          close_button = new QPushButton (this);
          close_button->setToolTip (tr ("Close ROI"));
          close_button->setIcon (QIcon (":/close.svg"));
          close_button->setEnabled (false);
          connect (close_button, SIGNAL (clicked()), this, SLOT (close_slot ()));
          layout->addWidget (close_button, 1);

          hide_all_button = new QPushButton (this);
          hide_all_button->setToolTip (tr ("Hide all ROIs"));
          hide_all_button->setIcon (QIcon (":/hide.svg"));
          hide_all_button->setCheckable (true);
          connect (hide_all_button, SIGNAL (clicked()), this, SLOT (hide_all_slot ()));
          layout->addWidget (hide_all_button, 1);

          main_box->addLayout (layout, 0);

          list_view = new QListView (this);
          list_view->setSelectionMode (QAbstractItemView::ExtendedSelection);
          list_view->setDragEnabled (true);
          list_view->viewport()->setAcceptDrops (true);
          list_view->setDropIndicatorShown (true);

          list_model = new ROI_Model (this);
          list_view->setModel (list_model);
          list_view->setSelectionMode (QAbstractItemView::SingleSelection);

          main_box->addWidget (list_view, 1);

          layout = new HBoxLayout;
          layout->setContentsMargins (0, 0, 0, 0);
          layout->setSpacing (0);

          draw_button = new QToolButton (this);
          draw_button->setToolButtonStyle (Qt::ToolButtonTextBesideIcon);
          QAction* action = new QAction (QIcon (":/draw.svg"), tr ("Edit"), this);
          action->setShortcut (tr ("E"));
          action->setToolTip (tr ("Add/remove voxels to/from ROI\n\nUse left mouse button to add voxels,\nright mouse button to erase"));
          action->setCheckable (true);
          action->setEnabled (false);
          connect (action, SIGNAL (toggled(bool)), this, SLOT (draw_slot ()));
          draw_button->setDefaultAction (action);
          layout->addWidget (draw_button, 1);

          undo_button = new QToolButton (this);
          undo_button->setToolButtonStyle (Qt::ToolButtonTextBesideIcon);
          action = new QAction (QIcon (":/undo.svg"), tr ("Undo"), this);
          action->setShortcut (tr ("Ctrl+Z"));
          action->setToolTip (tr ("Undo last edit"));
          action->setCheckable (false);
          action->setEnabled (false);
          connect (action, SIGNAL (triggered()), this, SLOT (undo_slot ()));
          undo_button->setDefaultAction (action);
          layout->addWidget (undo_button, 1);

          redo_button = new QToolButton (this);
          redo_button->setToolButtonStyle (Qt::ToolButtonTextBesideIcon);
          action = new QAction (QIcon (":/redo.svg"), tr ("Redo"), this);
          action->setShortcut (tr ("Ctrl+Y"));
          action->setToolTip (tr ("Redo last edit"));
          action->setCheckable (false);
          action->setEnabled (false);
          connect (action, SIGNAL (triggered()), this, SLOT (redo_slot ()));
          redo_button->setDefaultAction (action);
          layout->addWidget (redo_button, 1);

          main_box->addLayout (layout, 0);

          QGroupBox* group_box = new QGroupBox ("Draw mode");

          GridLayout* grid_layout = new GridLayout;
          group_box->setLayout (grid_layout);

          edit_mode_group = new QActionGroup (this);
          edit_mode_group->setExclusive (true);
          edit_mode_group->setEnabled (false);
          connect (edit_mode_group, SIGNAL (triggered (QAction*)), this, SLOT (select_edit_mode (QAction*)));

          rectangle_button = new QToolButton (this);
          rectangle_button->setToolButtonStyle (Qt::ToolButtonTextBesideIcon);
          action = new QAction (QIcon (":/rectangle.svg"), tr ("Rectangle"), this);
          action->setShortcut (tr ("Ctrl+R"));
          action->setToolTip (tr ("Edit ROI using a rectangle"));
          action->setCheckable (true);
          action->setChecked (false);
          edit_mode_group->addAction (action);
          rectangle_button->setDefaultAction (action);
          grid_layout->addWidget (rectangle_button, 0, 0);

          fill_button = new QToolButton (this);
          fill_button->setToolButtonStyle (Qt::ToolButtonTextBesideIcon);
          action = new QAction (QIcon (":/fill.svg"), tr ("Fill"), this);
          action->setShortcut (tr ("Ctrl+F"));
          action->setToolTip (tr ("Fill ROI slice"));
          action->setCheckable (true);
          action->setChecked (false);
          edit_mode_group->addAction (action);
          fill_button->setDefaultAction (action);
          grid_layout->addWidget (fill_button, 0, 1);

          brush_button = new QToolButton (this);
          brush_button->setToolButtonStyle (Qt::ToolButtonTextBesideIcon);
          action = new QAction (QIcon (":/brush.svg"), tr ("Brush"), this);
          action->setShortcut (tr ("Ctrl+B"));
          action->setToolTip (tr ("Edit ROI using a brush"));
          action->setCheckable (true);
          action->setChecked (true);
          edit_mode_group->addAction (action);
          brush_button->setDefaultAction (action);
          grid_layout->addWidget (brush_button, 0, 2);
          
          QLabel* label = new QLabel (tr("brush size: "));
          grid_layout->addWidget (label, 1, 0, 1, 2, Qt::AlignRight);

          brush_size_button = new AdjustButton (this);
          brush_size_button->setToolTip (tr ("Brush size (in mm)"));
          brush_size_button->setEnabled (true);
          grid_layout->addWidget (brush_size_button, 1, 2);

          main_box->addWidget (group_box, 0);

          layout = new HBoxLayout;
          layout->setContentsMargins (0, 0, 0, 0);
          layout->setSpacing (0);

          slice_copy_group = new QActionGroup (this);
          slice_copy_group->setEnabled (false);
          connect (slice_copy_group, SIGNAL (triggered (QAction*)), this, SLOT (slice_copy_slot (QAction*)));

          layout->addWidget (new QLabel ("Copy from slice: "), 0, Qt::AlignRight);

          copy_from_above_button = new QToolButton (this);
          copy_from_above_button->setToolButtonStyle (Qt::ToolButtonTextBesideIcon);
          action = new QAction (QIcon (":/copy_from_above.svg"), tr ("Above"), this);
          action->setToolTip (tr ("Copy data from the slice above into this slice"));
          action->setCheckable (false);
          action->setChecked (false);
          slice_copy_group->addAction (action);
          copy_from_above_button->setDefaultAction (action);
          layout->addWidget (copy_from_above_button, 1);

          copy_from_below_button = new QToolButton (this);
          copy_from_below_button->setToolButtonStyle (Qt::ToolButtonTextBesideIcon);
          action = new QAction (QIcon (":/copy_from_below.svg"), tr ("Below"), this);
          action->setToolTip (tr ("Copy data from the slice below into this slice"));
          action->setCheckable (false);
          action->setChecked (false);
          slice_copy_group->addAction (action);
          copy_from_below_button->setDefaultAction (action);
          layout->addWidget (copy_from_below_button, 1);

          main_box->addLayout (layout, 0);

          layout = new HBoxLayout;
          layout->setContentsMargins (0, 0, 0, 0);
          layout->setSpacing (0);

          colour_button = new QColorButton;
          colour_button->setEnabled (false);
          connect (colour_button, SIGNAL (clicked()), this, SLOT (colour_changed()));
          layout->addWidget (colour_button, 0);

          opacity_slider = new QSlider (Qt::Horizontal);
          opacity_slider->setToolTip (tr("ROI opacity"));
          opacity_slider->setRange (1,1000);
          opacity_slider->setSliderPosition (int (1000));
          connect (opacity_slider, SIGNAL (valueChanged (int)), this, SLOT (opacity_changed(int)));
          opacity_slider->setEnabled (false);
          layout->addWidget (opacity_slider, 1);

          main_box->addLayout (layout, 0);

          connect (list_view->selectionModel(),
              SIGNAL(selectionChanged(const QItemSelection &, const QItemSelection &)),
              SLOT (update_selection()));

          connect (&window(), SIGNAL (imageChanged()), this, SLOT (update_selection()));

          connect (list_model, SIGNAL (dataChanged (const QModelIndex&, const QModelIndex&)),
              this, SLOT (toggle_shown_slot (const QModelIndex&, const QModelIndex&)));

          update_selection();
        }






        ROI::~ROI()
        {
          for (int i = 0; i != list_model->rowCount(); ++i) {
            QModelIndex index = list_model->index (i, 0);
            ROI_Item* roi = list_model->get (index);
            if (!roi->saved) {
              if (QMessageBox::question (&window(), tr("ROI not saved"), tr (("Image " + roi->get_filename() + " has been modified. Do you want to save it?").c_str())) == QMessageBox::Yes)
                save (roi);
            }
          }
        }





        void ROI::new_slot ()
        {
          assert (window().image());
          list_model->create (window().image()->header());
          list_view->selectionModel()->clear();
          list_view->selectionModel()->select (list_model->index (list_model->rowCount()-1, 0, QModelIndex()), QItemSelectionModel::Select);
          updateGL ();
          in_insert_mode = false;
        }







        void ROI::open_slot ()
        {
          std::vector<std::string> names = Dialog::File::get_images (this, "Select ROI images to open");
          if (names.empty())
            return;
          std::vector<std::unique_ptr<MR::Image::Header>> list;
          for (size_t n = 0; n < names.size(); ++n)
            list.push_back (std::unique_ptr<MR::Image::Header> (new MR::Image::Header (names[n])));

          load (list);
          in_insert_mode = false;
        }





        void ROI::save (ROI_Item* roi)
        {
          std::vector<GLubyte> data (roi->info().dim(0) * roi->info().dim(1) * roi->info().dim(2));
          { 
            Window::GrabContext context; 
            roi->texture().bind();
            gl::PixelStorei (gl::PACK_ALIGNMENT, 1);
            gl::GetTexImage (gl::TEXTURE_3D, 0, GL_RED, GL_UNSIGNED_BYTE, (void*) (&data[0]));
          }

          try {
            MR::Image::Header header;
            header.info() = roi->info();
            header.set_ndim(3);
            header.datatype() = DataType::Bit;
            std::string name = GUI::Dialog::File::get_save_image_name (&window(), "Select name of ROI to save", roi->get_filename());
            if (name.size()) {
              MR::Image::Buffer<bool> buffer (name, header);
              roi->save (buffer.voxel(), data.data());
            }
          }
          catch (Exception& E) {
            E.display();
          }
          in_insert_mode = false;
        }





        int ROI::normal2axis (const Point<>& normal, const MR::Image::Transform& transform) const
        {
          float x_dot_n = std::abs (transform.image2scanner_dir (Point<> (1.0, 0.0, 0.0)).dot (normal));
          float y_dot_n = std::abs (transform.image2scanner_dir (Point<> (0.0, 1.0, 0.0)).dot (normal));
          float z_dot_n = std::abs (transform.image2scanner_dir (Point<> (0.0, 0.0, 1.0)).dot (normal));
          if (x_dot_n > y_dot_n)
            return x_dot_n > z_dot_n ? 0 : 2;
          else
            return y_dot_n > z_dot_n ? 1 : 2;
        }





        void ROI::save_slot ()
        {
          QModelIndexList indices = list_view->selectionModel()->selectedIndexes();
          assert (indices.size() == 1);
          ROI_Item* roi = dynamic_cast<ROI_Item*> (list_model->get (indices[0]));
          save (roi);
        }






        void ROI::load (std::vector<std::unique_ptr<MR::Image::Header>>& list) 
        {
          list_model->load (list);
          list_view->selectionModel()->clear();
          list_view->selectionModel()->select (list_model->index (list_model->rowCount()-1, 0, QModelIndex()), QItemSelectionModel::Select);
          updateGL ();
        }







        void ROI::close_slot ()
        {
          QModelIndexList indices = list_view->selectionModel()->selectedIndexes();
          assert (indices.size() == 1);
          ROI_Item* roi = dynamic_cast<ROI_Item*> (list_model->get (indices[0]));
          if (!roi->saved) {
            if (QMessageBox::question (&window(), tr("ROI not saved"), tr ("ROI has been modified. Do you want to save it?")) == QMessageBox::Yes)
              save_slot();
          }

          list_model->remove_item (indices.first());
          updateGL();
          in_insert_mode = false;
        }







        void ROI::draw_slot ()
        {
          if (draw_button->isChecked())
            grab_focus ();
          else
            release_focus ();
        }







        void ROI::undo_slot () 
        {
          QModelIndexList indices = list_view->selectionModel()->selectedIndexes();
          if (indices.size() != 1) {
            WARN ("FIXME: shouldn't be here!");
            return;
          }
          ROI_Item* roi = dynamic_cast<ROI_Item*> (list_model->get (indices[0]));

          roi->undo();
          update_undo_redo();
          updateGL();
          in_insert_mode = false;
        }







        void ROI::redo_slot () 
        {
          QModelIndexList indices = list_view->selectionModel()->selectedIndexes();
          if (indices.size() != 1) {
            WARN ("FIXME: shouldn't be here!");
            return;
          }
          ROI_Item* roi = dynamic_cast<ROI_Item*> (list_model->get (indices[0]));

          roi->redo();
          update_undo_redo();
          updateGL();
          in_insert_mode = false;
        }







        void ROI::slice_copy_slot (QAction* action)
        {
          QModelIndexList indices = list_view->selectionModel()->selectedIndexes();
          if (indices.size() != 1) {
            WARN ("FIXME: shouldn't be here!");
            return;
          }
          
          ROI_Item* roi = dynamic_cast<ROI_Item*> (list_model->get (indices[0]));

          const Projection* proj = window().get_current_mode()->get_current_projection();
          if (!proj) return;
          const Point<> current_origin = proj->screen_to_model (window().mouse_position(), window().focus());
          current_axis = normal2axis (proj->screen_normal(), roi->transform());
          current_slice = std::lround (roi->transform().scanner2voxel (current_origin)[current_axis]);

          roi->start (ROI_UndoEntry (*roi, current_axis, current_slice));

          const int source_slice = current_slice + ((action == copy_from_above_button->defaultAction()) ? 1 : -1);
          if (source_slice < 0 || source_slice >= roi->info().dim (current_axis))
            return;

          ROI_UndoEntry source (*roi, current_axis, source_slice);
          roi->current().copy (*roi, source);
          updateGL();
          in_insert_mode = false;
        }








        void ROI::select_edit_mode (QAction*) 
        {
          brush_size_button->setEnabled (brush_button->isChecked());
        }







        void ROI::hide_all_slot () 
        {
          updateGL();
          in_insert_mode = false;
        }






        void ROI::draw (const Projection& projection, bool is_3D, int, int)
        {
          if (is_3D) return;

          if (!is_3D) {
            // set up OpenGL environment:
            gl::Enable (gl::BLEND);
            gl::Disable (gl::DEPTH_TEST);
            gl::DepthMask (gl::FALSE_);
            gl::ColorMask (gl::TRUE_, gl::TRUE_, gl::TRUE_, gl::TRUE_);
            gl::BlendFunc (gl::SRC_ALPHA, gl::ONE_MINUS_SRC_ALPHA);
            gl::BlendEquation (gl::FUNC_ADD);
          }

          for (int i = 0; i < list_model->rowCount(); ++i) {
            if (list_model->items[i]->show && !hide_all_button->isChecked()) {
              ROI_Item* roi = dynamic_cast<ROI_Item*>(list_model->items[i].get());
              //if (is_3D) 
              //window.get_current_mode()->overlays_for_3D.push_back (image);
              //else
              roi->render (shader, projection, projection.depth_of (window().focus()));
            }
          }

          if (!is_3D) {
            // restore OpenGL environment:
            gl::Disable (gl::BLEND);
            gl::Enable (gl::DEPTH_TEST);
            gl::DepthMask (gl::TRUE_);
          }
        }







        void ROI::toggle_shown_slot (const QModelIndex& index, const QModelIndex& index2) 
        {
          if (index.row() == index2.row()) {
            list_view->setCurrentIndex(index);
          } 
          else {
            for (size_t i = 0; i < list_model->items.size(); ++i) {
              if (list_model->items[i]->show) {
                list_view->setCurrentIndex (list_model->index (i, 0));
                break;
              }
            }
          }
          updateGL();
          in_insert_mode = false;
        }






        void ROI::update_slot () 
        {
          updateGL();
        }






        void ROI::colour_changed () 
        {
          QModelIndexList indices = list_view->selectionModel()->selectedIndexes();
          for (int i = 0; i < indices.size(); ++i) {
            ROI_Item* roi = dynamic_cast<ROI_Item*> (list_model->get (indices[i]));
            QColor c = colour_button->color();
            roi->colour = { { GLubyte (c.red()), GLubyte (c.green()), GLubyte (c.blue()) } };
          }
          updateGL();
        }







        void ROI::opacity_changed (int)
        {
          QModelIndexList indices = list_view->selectionModel()->selectedIndexes();
          for (int i = 0; i < indices.size(); ++i) {
            ROI_Item* roi = dynamic_cast<ROI_Item*> (list_model->get (indices[i]));
            roi->alpha = opacity_slider->value() / 1.0e3f;
          }
          window().updateGL();
          in_insert_mode = false;
        }







        void ROI::update_undo_redo () 
        {
          QModelIndexList indices = list_view->selectionModel()->selectedIndexes();

          if (indices.size()) {
            ROI_Item* roi = dynamic_cast<ROI_Item*> (list_model->get (indices[0]));
            undo_button->defaultAction()->setEnabled (roi->has_undo());
            redo_button->defaultAction()->setEnabled (roi->has_redo());
          }
          else {
            undo_button->defaultAction()->setEnabled (false);
            redo_button->defaultAction()->setEnabled (false);
          }
        }






        void ROI::update_selection () 
        {
          if (!window().image()) {
            setEnabled (false);
            return;
          }
          else 
            setEnabled (true);

          QModelIndexList indices = list_view->selectionModel()->selectedIndexes();
          bool enable = window().image() && indices.size();

          opacity_slider->setEnabled (enable);
          save_button->setEnabled (enable);
          close_button->setEnabled (enable);
          draw_button->defaultAction()->setEnabled (enable);
          colour_button->setEnabled (enable);
          edit_mode_group->setEnabled (enable);
          slice_copy_group->setEnabled (enable);
          brush_size_button->setEnabled (enable && brush_button->isChecked());

          update_undo_redo();

          if (!indices.size()) {
            draw_button->defaultAction()->setChecked (false);
            return;
          }

          ROI_Item* roi = dynamic_cast<ROI_Item*> (list_model->get (indices[0]));
          colour_button->setColor (QColor (roi->colour[0], roi->colour[1], roi->colour[2]));
          opacity_slider->setValue (1.0e3f * roi->alpha);

          brush_size_button->setMin (roi->min_brush_size);
          brush_size_button->setMax (roi->max_brush_size);
          brush_size_button->setRate (0.1f * roi->min_brush_size);
          brush_size_button->setValue (roi->brush_size);
        }









        bool ROI::mouse_press_event () 
        { 
          if (in_insert_mode || window().modifiers() != Qt::NoModifier)
            return false;

          if (window().mouse_buttons() != Qt::LeftButton && window().mouse_buttons() != Qt::RightButton)
            return false;

          in_insert_mode = true;
          insert_mode_value = (window().mouse_buttons() == Qt::LeftButton);
          update_cursor();

          QModelIndexList indices = list_view->selectionModel()->selectedIndexes();
          if (indices.size() != 1) {
            WARN ("FIXME: shouldn't be here!");
            return false;
          }

          const Projection* proj = window().get_current_mode()->get_current_projection();
          if (!proj) 
            return false;
          current_origin =  proj->screen_to_model (window().mouse_position(), window().focus());
          window().set_focus (current_origin);
          prev_pos = current_origin;


          // figure out the closest ROI axis, and lock to it:
          ROI_Item* roi = dynamic_cast<ROI_Item*> (list_model->get (indices[0]));
          current_axis = normal2axis (proj->screen_normal(), roi->transform());

          // figure out current slice in ROI:
          current_slice = std::lround (roi->transform().scanner2voxel (current_origin)[current_axis]);

          // floating-point version of slice location to keep it consistent on
          // mouse move:
          Point<> slice_axis (0.0, 0.0, 0.0);
          slice_axis[current_axis] = current_axis == 2 ? 1.0 : -1.0;
          slice_axis = roi->transform().image2scanner_dir (slice_axis);
          current_slice_loc = current_origin.dot (slice_axis);

          Math::Versor<float> orient;
          orient.from_matrix (roi->info().transform());
          window().set_snap_to_image (false);
          window().set_orientation (orient);

          roi->start (ROI_UndoEntry (*roi, current_axis, current_slice));
         

          if (brush_button->isChecked()) {
            if (brush_size_button->isMin())
              roi->current().draw_line (*roi, prev_pos, current_origin, insert_mode_value);
            else
              roi->current().draw_circle (*roi, current_origin, insert_mode_value, brush_size_button->value());
          } else if (rectangle_button->isChecked()) {
            roi->current().draw_rectangle (*roi, current_origin, current_origin, insert_mode_value);
          } else if (fill_button->isChecked()) {
            roi->current().draw_fill (*roi, current_origin, insert_mode_value);
          }


          updateGL();

          return true; 
        }






        bool ROI::mouse_move_event () 
        { 
          if (!in_insert_mode)
            return false;

          QModelIndexList indices = list_view->selectionModel()->selectedIndexes();
          if (!indices.size()) {
            WARN ("FIXME: shouldn't be here!");
            return false;
          }
          ROI_Item* roi = dynamic_cast<ROI_Item*> (list_model->get (indices[0]));

          const Projection* proj = window().get_current_mode()->get_current_projection();
          if (!proj) 
            return false;

          Point<> pos =  proj->screen_to_model (window().mouse_position(), window().focus());
          Point<> slice_axis (0.0, 0.0, 0.0);
          slice_axis[current_axis] = current_axis == 2 ? 1.0 : -1.0;
          slice_axis = roi->transform().image2scanner_dir (slice_axis);
          float l = (current_slice_loc - pos.dot (slice_axis)) / proj->screen_normal().dot (slice_axis);
          window().set_focus (window().focus() + l * proj->screen_normal());
          const Point<> pos_adj = pos + l * proj->screen_normal();

          if (brush_button->isChecked()) {
            if (brush_size_button->isMin())
              roi->current().draw_line (*roi, prev_pos, pos_adj, insert_mode_value);
            else {
              const float diameter = brush_size_button->value();
              roi->current().draw_thick_line (*roi, prev_pos, pos_adj, insert_mode_value, diameter);
              roi->current().draw_circle (*roi, pos_adj, insert_mode_value, diameter);
            }
          } else if (rectangle_button->isChecked()) {
            roi->current().draw_rectangle (*roi, current_origin, pos_adj, insert_mode_value);
          } else if (fill_button->isChecked()) {
            // Do nothing
          }

          updateGL();
          prev_pos = pos_adj;
          return true; 
        }






        bool ROI::mouse_release_event () 
        { 
          in_insert_mode = false;
          update_cursor();
          update_undo_redo();
          return true; 
        }






        QCursor* ROI::get_cursor ()
        {
          if (!draw_button->isChecked())
            return nullptr;
          if (in_insert_mode && !insert_mode_value)
            return &Cursor::erase;
          return &Cursor::draw;
        }







        void ROI::add_commandline_options (MR::App::OptionList& options) 
        { 
          using namespace MR::App;
          options
            + OptionGroup ("ROI editor tool options")

            + Option ("roi.load", "Loads the specified image on the ROI editor tool.")
            +   Argument ("image").type_image_in()

            + Option ("roi.opacity", "Sets the overlay opacity to floating value [0-1].")
            +   Argument ("value").type_float (0.0, 1.0, 1.0);
        }






        bool ROI::process_commandline_option (const MR::App::ParsedOption& opt) 
        {
          if (opt.opt->is ("roi.load")) {
            std::vector<std::unique_ptr<MR::Image::Header>> list;
            try { list.push_back (std::unique_ptr<MR::Image::Header> (new MR::Image::Header (opt[0]))); }
            catch (Exception& e) { e.display(); }
            load (list);
            return true;
          }

          if (opt.opt->is ("roi.opacity")) {
            try {
              float value = opt[0];
              opacity_slider->setSliderPosition(int(1.e3f*value));
            }
            catch (Exception& e) { e.display(); }
            return true;
          }

          return false;
        }



      }
    }
  }
}




