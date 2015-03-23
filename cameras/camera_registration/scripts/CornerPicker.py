#!/usr/bin/env python

import sys
import os
import pickle
import math

import wx
import cv2
import numpy as np

class AppFrame(wx.Frame):

    def __init__(self, parent, img, corners, pattern_size):

        super(AppFrame, self).__init__(parent,
                          title="CornerPicker",
                          size=(800,650))


        self._orig_img = img
        self._orig_corners = corners
        self.corners = corners.copy()
        self.pattern_size = pattern_size
        if not corners.shape[0] == (pattern_size[0] * pattern_size[1]):
            print "Corner size doesn't match pattern size!"
            print "Reseting corners to defaults"
            self._reset_corners()
        self._init_ui()
        self._display_image()

    def _init_ui(self):

        self._setup_toolbar()
        self._setup_panels()

        self.Center()
        self.Show()

    def _setup_toolbar(self):

        self.toolbar = self.CreateToolBar()

        art_bmp = wx.ArtProvider.GetBitmap
        open_tool = self.toolbar.AddLabelTool(wx.ID_OPEN, "Open",
                art_bmp(wx.ART_FILE_OPEN, wx.ART_TOOLBAR, (15,15)),
                shortHelp="Open Boards")
        save_tool = self.toolbar.AddLabelTool(wx.ID_SAVE, "Save",
                    art_bmp(wx.ART_FILE_SAVE, wx.ART_TOOLBAR, (15, 15)),
                    shortHelp="Save Coordinates")

        self.toolbar.Realize()

#        self.Bind(wx.EVT_TOOL, self._open, open_tool)
        #self.Bind(wx.EVT_TOOL, self._save_coordinates, save_tool)

    def _setup_panels(self):

        vbox = wx.BoxSizer(wx.VERTICAL)
        self.impanel = wx.Panel(self, size=(640, 480), style=wx.NO_BORDER)
        self.impanel.SetBackgroundColour('WHITE')
        self.bmp_display = wx.StaticBitmap(self.impanel, size=(640,480))
        self.image_text = wx.StaticText(self, wx.ID_ANY, "Image Name")

        vbox.Add(self.impanel, 0, wx.TOP | wx.ALIGN_CENTER, border=10)

        self.button_panel = wx.Panel(self, size=(700, 60), style=wx.RAISED_BORDER)
        self.button_panel.SetBackgroundColour('#53728c')

        but_vbox = wx.BoxSizer(wx.VERTICAL)

        but_hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        but_hbox2 = wx.BoxSizer(wx.HORIZONTAL)

        next_button = wx.Button(self.button_panel, label="Next")
        prev_button = wx.Button(self.button_panel, label="Previous")

        but_hbox2.Add(prev_button, 0, wx.ALIGN_CENTER | wx.LEFT| wx.RIGHT, border=10)
        but_hbox2.Add(next_button, 0, wx.ALIGN_CENTER | wx.LEFT| wx.RIGHT, border=10)

        but_vbox.Add(but_hbox1, 0, wx.ALIGN_CENTER)
        but_vbox.Add(but_hbox2, 0, wx.ALIGN_CENTER)

        self.button_panel.SetSizer(but_vbox)

        vbox.Add(self.image_text, 0, wx.ALL | wx.ALIGN_CENTER, border = 10)
        vbox.Add(self.button_panel, 0, wx.ALL | wx.ALIGN_CENTER, border=10)

        self.SetSizer(vbox)

        # Bind button events to callbacks
#        self.Bind(wx.EVT_BUTTON, self._next_callback, next_button)
#        self.Bind(wx.EVT_BUTTON, self._prev_callback, prev_button)
        self.bmp_display.Bind(wx.EVT_LEFT_DOWN, self._on_mouse_button)
        self.Bind(wx.EVT_CLOSE, self._on_close)

    ####
    #def _open(self, event):
    #    #fd = wx.FileDialog(None, "Choose the image to open", "../..", style=wx.FD_OPEN)
    #    #if fd.ShowModal() == wx.ID_OK:
    #    full_image_path = fd.GetPath()
    #    (path, image_name) = os.path.split(full_image_path)
    #    os.chdir(path)
    #    self.current_image_name = image_name
    #    self.E_coords = ECoords(os.path.join(path, "E_coords.pydict"))
    #    #else:
    #    #    full_image_path = None
    #    #    image_name = None

    #    #fd.Destroy()
    #    #import pdb; pdb.set_trace()
    #    #print "You pressed open!", image_name
    #    #self._create_image_lists()
    #    self._display_image(image_name)

    #def _save_coordinates(self, event):
    #    fd = wx.FileDialog(None, "Save the current coordinates...", ".", "E_coords.pydict", style=wx.FD_SAVE)
    #    if fd.ShowModal() == wx.ID_OK:
    #        file_path = fd.GetPath()
    #        f = open(file_path, 'w')
    #        pickle.dump(self.E_coords.coords, f)
    #        self.E_coords.is_saved = True
    #        print "You saved the coordinates"
    #    else:
    #        print "Save canceled"

    #    fd.Destroy()
    #"""
    def _display_image(self):

        # Reset img to a copy of the original
        img = self._draw_corners_on_image()
        image_buffer = img.tostring()
        wximg = wx.ImageFromData(640, 480, image_buffer)
        bmp = wx.BitmapFromImage(wximg)
        self.bmp_display.SetBitmap(bmp)
        image_name = "Unfound Corners"
        self.current_image_name = image_name
        self.image_text.SetLabel(image_name)
        print "Displaying Image %s!" % image_name
# TODO Gonna have to update the image with new drawn points once we click on it...

    def _on_mouse_button(self, event):
        """
        Move the nearest point to the cursor position.
        """
        pt = event.GetPositionTuple()
        print pt
        self._replace_nearest_point(pt)
        self._display_image()


    def _on_right_mouse_button(self, event):
        """
        Remove the nearest point and all following points
        """
        pass

    def _replace_nearest_point(self, click_point):
        """
        Find the point nearest to the click and replace it with the click value
        """
		
		
        # Find the nearest point inside the corners matrix
        min_distance = 10000
        dist_index = None
        for corner_point, index in zip(self.corners, xrange(self.corners.shape[0])):
            distance = math.sqrt((corner_point[0] - click_point[0]) **2 +
                                (corner_point[1] - click_point[1]) **2)
            if distance < min_distance:
                min_distance = distance
                dist_index = index

        self.corners[dist_index] = list(click_point)

        f = open('corners_tmp.pkl', 'w')
        pickle.dump(self.corners, f)
        f.close()


    def _draw_corners_on_image(self):

        img = self._orig_img.copy()
        radius = 7
        color = (255,0,0)
        # Draw lines between adjacent points to visualize correct order
        for pt1, pt2 in zip(self.corners[:-1], self.corners[1:]):
            cv2.line(img, tuple(pt1), tuple(pt2), color)

        for pt in self.corners:
            cv2.circle(img, tuple(pt), radius, color)
            start1 = (pt[0] - radius, pt[1])
            end1 = (pt[0]+radius, pt[1])
            start2 = (pt[0], pt[1]-radius)
            end2 = (pt[0], pt[1]+radius)
            cv2.line(img, start1, end1, color)
            cv2.line(img, start2, end2, color)

        return img

    def _reset_corners(self):
        num_rows = self.pattern_size[0]
        num_cols = self.pattern_size[1]
        x_pts = np.r_[[[x]*num_rows for x in range(11,0,-1)]].flatten()
        y_pts = np.array(range(1,5)*num_cols)
        corners = np.c_[x_pts, y_pts]
        corners *= 40
        corners = np.array(corners, dtype='float32')
        self.corners = corners

# TODO Fill out the mouse event and add a right click one to remove points

    def _on_close(self, event):
        #if not self.E_coords.is_saved:
         #   self._save_coordinates(event)

        if event.CanVeto():
            self.Destroy()
        else:
            self.Destroy()
    #"""

# TODO Make this appropriate for the corners
class ECoords(object):

    def __init__(self, path_to_data=''):

        if os.path.isfile(path_to_data):
            with open(path_to_data, 'r') as f:
                self._coord_data = pickle.load(f)
        else:
            self._coord_data = dict()

        self.is_saved = True

    @property
    def coords(self):
        return self._coord_data

    @coords.setter
    def coords(self, key, value):
        print "Setting"
        self._coord_data[key] = value
        self.is_saved = False

def main(corner_image, corners, pattern_size):
    print corners
    print corners.shape, pattern_size
    f = open('corners_tmp.pkl', 'w')
    pickle.dump(corners, f)
    f.close()
    app = wx.App()
    main_frame = AppFrame(None, corner_image, corners, pattern_size)
    app.MainLoop()
    f = open('corners_tmp.pkl', 'r')
    corners = pickle.load(f)
    f.close()
    os.remove('corners_tmp.pkl')
    print corners
    return corners

if __name__=="__main__":
    main()
