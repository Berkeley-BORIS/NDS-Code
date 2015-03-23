#!/usr/bin/env python

import sys
import os
import pickle
import fnmatch

import wx
import cv2

class AppFrame(wx.Frame):

    def __init__(self, parent):

        super(AppFrame, self).__init__(parent,
                          title="EPicker",
                          size=(800,650))

        self._init_ui()

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

        self.Bind(wx.EVT_TOOL, self._open, open_tool)
        self.Bind(wx.EVT_TOOL, self._save_coordinates, save_tool)

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
        self.Bind(wx.EVT_BUTTON, self._next_callback, next_button)
        self.Bind(wx.EVT_BUTTON, self._prev_callback, prev_button)
        self.bmp_display.Bind(wx.EVT_LEFT_DOWN, self._on_left_mouse_button)
        self.bmp_display.Bind(wx.EVT_RIGHT_DOWN, self._on_right_mouse_button)
        self.Bind(wx.EVT_CLOSE, self._on_close)

    def _open(self, event):
        default_path = "../../image_rectification/data/"  # "../../"
        #default_path = "/Users/natdispstats/Documents/cameras/disparity_estimation/data/"
        fd = wx.FileDialog(None, "Choose the image to open", default_path, style=wx.FD_OPEN)
        if fd.ShowModal() == wx.ID_OK:
            full_image_path = fd.GetPath()
            (path, image_name) = os.path.split(full_image_path)
            print path
            if fnmatch.fnmatch(path, '*_rect'):
                self.is_rectified = True
            else:
                self.is_rectified = False
            os.chdir(path)
            self.current_image_name = image_name
            self.E_coords = ECoords(os.path.join(path, "E_coords.pydict"))
        else:
            full_image_path = None
            image_name = None

        fd.Destroy()
        print "You pressed open!", image_name
        if fnmatch.fnmatch(image_name, '*.bmp'):
            if self.get_frames_file(path):
                self._create_image_lists_from_data()
            else:
                self._create_image_lists()
        elif fnmatch.fnmatch(image_name, '*_frames.txt'):
            self.frames_fpath = os.path.join(path, image_name)
            self._create_image_lists_from_data()
        else:
            raise ValueError("Unknown how to process that file!")
        self._display_image(self.current_image_name)

    def get_frames_file(self, path):
        """
        Searches in current directory and a few other places if they exist for
        a frames.txt file
        """

        files = os.listdir('.')
        frames_file = fnmatch.filter(files, '*_frames.txt')
        if len(frames_file) > 0:
            self.frames_fpath = frames_file[0]
            return True

        root, path_dname = os.path.split(path)
        split_path_dname = path_dname.split("_")
        if split_path_dname[0] == "radials":
            root, img_dname = os.path.split(root)
            trial_ID = img_dname.split("_")[0]
        else:
            trial_ID = split_path_dname[0]
            frames_fname = "_".join(split_path_dname[:-1]) + ".txt"

        frames_fpath = os.path.join("/Users", "natdispstats", "Documents", "cameras", "e_picking", "data", trial_ID, frames_fname)
        print frames_fpath
        if os.path.isfile(frames_fpath):
            self.frames_fpath = frames_fpath
            return True

        return False


    def _save_coordinates(self, event):
        fd = wx.FileDialog(None, "Save the current coordinates...", ".", "E_coords.pydict", style=wx.FD_SAVE)
        if fd.ShowModal() == wx.ID_OK:
            file_path = fd.GetPath()
            f = open(file_path, 'w')
            pickle.dump(self.E_coords.coords, f)
            self.E_coords.is_saved = True
            print "You saved the coordinates"
        else:
            print "Save canceled"

        fd.Destroy()

    def _display_image(self, image_name):

        #loaded_im = cv2.imread(image_name,flags=1)

        bmp = wx.Bitmap(image_name, type=wx.BITMAP_TYPE_BMP)

        print "Displaying Image %s!" % image_name
        components = image_name[:-4].split('_')
        if self.is_rectified:
            frame_num = int(components[3])  # NOTE: This assumes image names start with "rect"
        else:
            frame_num = int(components[2])

        if not self.E_coords.coords.get(frame_num) == None:
            value = self.E_coords.coords[frame_num]
            pen = wx.Pen('blue',1)
            dc = wx.MemoryDC()
            dc.SelectObject(bmp)
            #dc.Clear()
            dc.SetPen(pen)
            dc.DrawCircle(int(value[0]),int(value[1]),2)
            dc.SelectObject(wx.NullBitmap)
            self.Show(True)
            #self.Drawcircle(value[0],value[1],10)

        self.bmp_display.SetBitmap(bmp)
        self.current_image_name = image_name
        self.image_text.SetLabel(image_name)

        try:
            value = self.E_coords.coords[frame_num]
            print frame_num, value
            self.SetBackgroundColour(wx.RED)

        except KeyError:
            self.SetBackgroundColour(wx.LIGHT_GREY)

    def _create_image_lists(self):

        image_list = os.listdir('.')
        self._cam1_images = {}
        self._cam2_images = {}
        self._cam1_images['smallest_frame_num'] = '9999999'
        self._cam1_images['largest_frame_num'] = '0'
        self._cam2_images['smallest_frame_num'] = '9999999'
        self._cam2_images['largest_frame_num'] = '0'
        if self.is_rectified:
            cam_index = 1
            frame_index = 3
        else:
            cam_index = 0
            frame_index = 2
        for image in image_list:
            if image.find('cam') > -1:
                components = image[:-4].split("_")
                frame_num = components[frame_index]
                if components[cam_index] == "cam1":
                    self._cam1_images[frame_num] = image
                    self._cam1_images[frame_num] = image
                elif components[cam_index] == "cam2":
                    self._cam2_images[frame_num] = image

                # update smallest and largest frame_nums for cam1
                if int(frame_num) > int(self._cam1_images['largest_frame_num']):
                    self._cam1_images['largest_frame_num'] = frame_num
                if int(frame_num) < int(self._cam1_images['smallest_frame_num']):
                    self._cam1_images['smallest_frame_num'] = frame_num

                # update smallest and largest frame_nums for cam2
                if int(frame_num) > int(self._cam2_images['largest_frame_num']):
                    self._cam2_images['largest_frame_num'] = frame_num
                if int(frame_num) < int(self._cam2_images['smallest_frame_num']):
                    self._cam2_images['smallest_frame_num'] = frame_num

     #   for image in image_list:
     #       if image.find('cam') > -1:
     #           components = image[:-4].split("_")
     #           if components[cam_index] == "cam1":
     #               frame_index = int(components[3])
     #               self._cam1_images[frame_index] = image
     #           elif components[cam_index] == "cam2":
     #               frame_index = int(components[3])
     #               self._cam2_images[frame_index] = image

    def _create_image_lists_from_data(self):

        print "Using", self.frames_fpath, "to filter frame list!"
        self._cam1_images = {}
        self._cam1_images['smallest_frame_num'] = '9999999'
        self._cam1_images['largest_frame_num'] = '0'
        with open(self.frames_fpath, 'r') as frame_file:
            for frame_num in frame_file:
                frame_num = int(frame_num.strip())
                if self.is_rectified:
                    frame_num = "%06d" % frame_num
                    image = "_".join(['rect', 'cam1', 'frame', frame_num]) + ".bmp"
                else:
                    frame_num = str(frame_num)
                    image = "_".join(['cam1', 'frame', frame_num]) + ".bmp"
                self._cam1_images[frame_num] = image

                # update smallest and largest frame_nums
                if int(frame_num) > int(self._cam1_images['largest_frame_num']):
                    self._cam1_images['largest_frame_num'] = frame_num
                if int(frame_num) < int(self._cam1_images['smallest_frame_num']):
                    self._cam1_images['smallest_frame_num'] = frame_num
        self.current_image_name = self._cam1_images[self._cam1_images['smallest_frame_num']]

    def _prev_callback(self, event):
        image_name = self._find_image(-1)
        self._display_image(image_name)
        print "Previous!"

    def _next_callback(self, event):
        image_name = self._find_image(1)
        self._display_image(image_name)
        print "Next!"

    def _find_image(self, direction):

        current_image_name = self.current_image_name
        components = current_image_name[:-4].split("_")
        if self.is_rectified:
            cam_index = 1
            frame_index = 3
        else:
            cam_index = 0
            frame_index = 2

        if components[cam_index] == "cam1":
            list_to_use = self._cam1_images
        elif components[cam_index] == "cam2":
            list_to_use = self._cam2_images
        frame = int(components[frame_index])

        while True:
                next_frame = frame + direction
                if self.is_rectified:
                    next_frame = "%06d" % next_frame
                else:
                    next_frame = str(next_frame)

                if list_to_use.has_key(next_frame):
                    next_image_name = list_to_use[next_frame]
                    break

                if int(next_frame) > int(list_to_use['largest_frame_num']):
                    next_image_name = list_to_use[list_to_use['smallest_frame_num']]
                    break
                elif int(next_frame) < int(list_to_use['smallest_frame_num']):
                    next_image_name = list_to_use[list_to_use['largest_frame_num']]
                    break

                frame += direction

        return next_image_name

    def _on_left_mouse_button(self, event):
        pt = event.GetPositionTuple()
        print pt
        components = self.current_image_name[:-4].split("_")
        if self.is_rectified:
            frame_index = 3
        else:
            frame_index = 2
        frame_num = int(components[frame_index])
        self.E_coords.coords[frame_num] = pt
        self.E_coords.is_saved = False
        print "Set", frame_num, "to", pt
        self._next_callback(event)

    def _on_right_mouse_button(self, event):
        components = self.current_image_name[:-4].split("_")
        if self.is_rectified:
            frame_index = 3
        else:
            frame_index = 2
        frame_num = int(components[frame_index])
        if self.E_coords.coords.has_key(frame_num):
            val = self.E_coords.coords.pop(frame_num)
            print "Frame", frame_num, "removed! Had value", val
            self.E_coords.is_saved = False
        else:
            print "No frame", frame_num, "to remove!"
        self.SetBackgroundColour(wx.LIGHT_GREY)
        self._display_image(self.current_image_name)

    def _on_key_press(self, event):
        print "A key!", event.GetKeyCode()
        if event.GetKeyCode() == wx.WXK_LEFT:
            self._prev_callback(event)
        elif event.GetKeyCode() == wx.WXK_RIGHT:
            self._next_callback(event)
        event.Skip()

    def _on_close(self, event):
        if not self.E_coords.is_saved:
            self._save_coordinates(event)

        if event.CanVeto():
            self.Destroy()
        else:
            self.Destroy()

class ECoords(object):

    def __init__(self, path_to_data=''):

        if os.path.isfile(path_to_data):
            with open(path_to_data, 'r') as f:
                self._coord_data = pickle.load(f)
        else:
            self._coord_data = dict()

        self.is_saved = True
        self.background = wx.BLACK
        self.coords = self._coord_data

    # TODO These properties are not working! Watch out!
    """
    def get_coords(self):
        print "Getter invoked"
        return self._coord_data

    def set_coords(self, key, value):

        print "Setter invoked"
        self.background = wx.BLUE
        # Remove the key entirely if the value is "remove"
        if value == 'remove':
            self.background = wx.BLUE
            if self._coord_data.has_key(key):
                orig_value = self._coord_data.pop(key)
                print "Removing value", orig_value, "with key", key
            else:
                print "No key", key, "to pop!"
        else:  # Set the value as normal
            print "Setting"
            self._coord_data[key] = value

        self.is_saved = False
    """

def main():
    args = sys.argv
    app = wx.App()
    main_frame = AppFrame(None)
    app.MainLoop()

if __name__=="__main__":
    main()
