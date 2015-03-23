import sys
import os.path

import pygame
import cv2
import cv
import scipy as sp

def play_movie(image_dir, fixation_file=None, fps=30):

    pygame.init()

    clock = pygame.time.Clock()
    flags = pygame.NOFRAME
    depth = 32
    surf = pygame.display.set_mode((640,480), flags, depth)

    im_base_name = "cam1_frame_"
    im_extension = ".bmp"

    if fixation_file is not None:
        fixations = sp.load(fixation_file)
        fixations[sp.isnan(fixations)] = -100
        fixations[abs(fixations) > 1000] = 1000
    else:
        fixations = []

    try:
        pygame.event.clear()
        pygame.event.set_allowed(None)
        pygame.event.set_allowed(pygame.KEYDOWN)
        for framenum in xrange(4000):
            im_name = "".join([im_base_name, str(framenum), im_extension])
            im_path = os.path.join(image_dir, im_name)
            im = cv2.imread(im_path)
            if len(fixations) != 0:
                if sp.floor(fixations[framenum, 1]) == 237:
                    continue
                cv2.circle(im, tuple(fixations[framenum]), 3, (255, 255, 255))
            im_buf = im.tostring()
            im = pygame.image.frombuffer(im_buf, (640,480), "RGB")
            surf.blit(im, (0,0))
            pygame.display.flip()
            print "Frame", framenum
            for event in pygame.event.get():
                if event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_ESCAPE:
                        sys.exit(0)

            clock.tick(30)
    except KeyboardInterrupt:
        print "Quitting"
    except Exception as e:
        print "An exception!"
        print e
        raise
    finally:
        pygame.quit()

def make_video(image_dir, filename="vidout.avi", fixation_file=None):

    MPEG_FOURCC = 827148624
    vwriter = cv2.VideoWriter()

    if fixation_file is not None:
        fixations = sp.load(fixation_file)
        fixations[sp.isnan(fixations)] = -100
        fixations[abs(fixations) > 1000] = 1000
    else:
        fixations = []

    im_base_name = "cam1_frame_"
    im_extension = ".bmp"

    suc = vwriter.open(os.path.join(image_dir, filename), cv.CV_FOURCC('M', 'J', 'P', 'G'), 30, (640,480))

    if not suc:
        raise IOError("Failed to open movie")

    for frame_num in xrange(1000):
        im_name = "".join([im_base_name, str(frame_num), im_extension])
        im_path = os.path.join(image_dir, im_name)
        im = cv2.imread(im_path)
        if len(fixations) != 0:
            cv2.circle(im, tuple(fixations[frame_num]), 3, (255, 255, 255))
        vwriter.write(im)


if __name__=="__main__":
    image_dir = sys.argv[1]
    try:
        fixation_file = sys.argv[2]
    except IndexError:
        fixation_file = None

    play_movie(image_dir, fixation_file)
