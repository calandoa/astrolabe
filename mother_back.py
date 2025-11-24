#!/usr/bin/python3
# mother_back.py
# -*- coding: utf-8 -*-
#
# The python script in this file makes the various parts of a model astrolabe.
#
# Copyright (C) 2010-2024 Dominic Ford <https://dcford.org.uk/>
#
# This code is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# You should have received a copy of the GNU General Public License along with
# this file; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA  02110-1301, USA

# ----------------------------------------------------------------------------

"""
Render the back of the mother of the astrolabe.
"""

from math import pi, sin, tan, cos, acos, atan, atan2, floor
from typing import Dict

from constants import unit_deg, unit_rev, unit_cm, unit_mm, centre_scaling, r_1, d_12
from graphics_context import BaseComponent
from numpy import arange
from settings import fetch_command_line_arguments
from text import text
from themes import themes
import calendar
#import scipy.interpolate
from datetime import datetime, timedelta


class MotherBack(BaseComponent):
    """
    Render the back of the mother of the astrolabe.
    """

    def default_filename(self) -> str:
        """
        Return the default filename to use when saving this component.
        """
        return "mother_back"

    def bounding_box(self, settings: dict) -> Dict[str, float]:
        """
        Return the bounding box of the canvas area used by this component.

        :param settings:
            A dictionary of settings required by the renderer.
        :return:
            Dictionary with the elements 'x_min', 'x_max', 'y_min' and 'y_max' set
        """

        r_outer: float = r_1 + 0.4 * unit_cm

        return {
            'x_min': -r_outer,
            'x_max': r_outer,
            'y_min': -r_outer - 2 * unit_cm,
            'y_max': r_outer
        }

    def do_rendering(self, settings, context):
        """
        This method is required to actually render this item.

        :param settings:
            A dictionary of settings required by the renderer.
        :param context:
            A GraphicsContext object to use for drawing
        :return:
            None
        """

        language = settings['language']
        theme = themes[settings['theme']]

        context.set_color(color=theme['lines'])

        # Define the radii of all the concentric circles to be drawn on back of mother

        d_12 = 0.03 * r_1
        dr = d_12 * 3 / 4

        # Scale of angles around the rim of the astrolabe
        r_2 = r_1 - d_12
        r_3 = r_2 - d_12

        # Zodiacal constellations
        #r_4 = r_3 - d_12
        r_5 = r_3 - d_12*1.2

        # Calendar for 1394
        #r_6 = r_5 - d_12
        #r_7 = r_6 - d_12

        # Days of the year
        #r_8 = r_5 - d_12 / 2

        # Calendar for 1974
        r_9 = r_5 - 4*dr
        r_10 = r_9 - d_12

        # Saints' days
        #r_11 = r_10 - d_12
        #r_12 = r_11 - d_12
        r_12 = r_5 - d_12 * 5

        # Radius of the central hole
        r_13 = d_12 * centre_scaling

        ## Draw the handle at the top of the astrolabe
        #ang = 180 * unit_deg - acos(unit_cm / r_1)
        #context.begin_path()
        #context.arc(centre_x=0, centre_y=-r_1, radius=2 * unit_cm,
        #            arc_from=-ang - pi / 2, arc_to=ang - pi / 2)
        ##context.move_to(x=0, y=-r_1 - 2 * unit_cm)
        ##context.line_to(x=0, y=-r_1 + 2 * unit_cm)
        #context.stroke()

        # Draw circles 1-13 onto back of mother
        context.begin_path()
        context.circle(centre_x=0, centre_y=0, radius=r_1)
        context.begin_sub_path()
        context.circle(centre_x=0, centre_y=0, radius=r_13)
        context.stroke(line_width=1)
        context.clip()

        for radius, line_width in ((r_2, 2), (r_3, 1), (r_5, 2), #(r_8, 1), (r_9, 1),  # (r_6, 1), (r_7, 1),
                                   (r_12, 1), (r_13, 1)): # (r_10, 2), (r_11, 1),
            context.begin_path()
            context.circle(centre_x=0, centre_y=0, radius=radius)
            context.stroke(line_width=line_width)

        # Label space between circles 1-5 with passage of Sun through zodiacal constellations

        # Mark every 30 degrees, where Sun enters new zodiacal constellation
        #for theta in arange(0 * unit_deg, 359 * unit_deg, 30 * unit_deg):
        #    context.begin_path()
        #    context.move_to(x=r_2 * cos(theta), y=-r_2 * sin(theta))
        #    context.line_to(x=r_5 * cos(theta), y=-r_5 * sin(theta))
        #    context.stroke(line_width=2)

        ## Mark 5-degree intervals within each zodiacal constellation
        #for theta in arange(0 * unit_deg, 359 * unit_deg, 5 * unit_deg):
        #    context.begin_path()
        #    context.move_to(x=r_2 * cos(theta), y=-r_2 * sin(theta))
        #    context.line_to(x=r_4 * cos(theta), y=-r_4 * sin(theta))
        #    context.stroke(line_width=1)

        # Mark fine scale of 1-degree intervals between circles 2 and 3
        for theta in range(360):
            if theta % 30 == 0:
                line_width = 2
                rv = r_5
            elif theta % 10 == 0:
                line_width = 1
                rv = r_3
            elif theta % 5 == 0:
                line_width = 1
                rv = (r_2 + 3*r_3)/4
            else:
                line_width = 1
                rv = (r_2 + r_3)/2

            context.begin_path()
            context.move_to(x=r_2 * cos(theta*unit_deg), y=-r_2 * sin(theta*unit_deg))
            context.line_to(x=rv * cos(theta*unit_deg), y=-rv * sin(theta*unit_deg))
            context.stroke(line_width = line_width)

        # Between circles 1 and 2, surround the entire astrolabe with a protractor scale from 0 to 90 degrees

        # Radius from centre for writing the text of the protractor scale around the rim
        rt_1 = (r_1 + r_2) / 2 - d_12 * .15

        for theta in arange(-180 * unit_deg, 179 * unit_deg, 10 * unit_deg):
            # Work out angle to display around the rim: counts from 0 to 90 four times, not -180 to 180 degrees!
            if theta < -179 * unit_deg:
                theta_disp = 0
            elif theta < - 90 * unit_deg:
                theta_disp = theta + 180 * unit_deg
            elif theta < 0 * unit_deg:
                theta_disp = -theta
            elif theta < 90 * unit_deg:
                theta_disp = theta
            else:
                theta_disp = -theta + 180 * unit_deg

            # Display angles around rim as rounded integers
            theta_disp = floor(theta_disp / unit_deg + 0.01)

            context.set_font_size(.8)

            theta2 = theta #+ 0.2 * unit_deg
            context.text(text="{:.0f}".format(theta_disp),
                         x=rt_1 * cos(theta2), y=-rt_1 * sin(theta2),
                         h_align=0, v_align=.4, gap=0, rotation=-theta - 90 * unit_deg)

        # Between circles 3 and 4, mark 10-, 20-, 30-degree points within each zodiacal constellation

        # Radius for writing the 30 degree scales within each zodiacal constellation
        rt_2 = (r_3 * 0.8 + r_5 * 0.2)
        context.set_font_size(.6)

        for theta in range(-175, 176, 10):
            theta_disp = (theta + 360) % 30
            if theta_disp == 15:
                continue
            # Work out what angle to display, which is rotation angle modulo 30 degrees

            # Write two digits separately, with a slight gap between them for the dividing line they label
            context.text(text="{:.0f}".format(theta_disp),
                         x=rt_2 * cos(theta * unit_deg), y=-rt_2 * sin(theta* unit_deg),
                         h_align=0, v_align=-0.8, gap=0, rotation= -(theta+90)* unit_deg )

        context.set_font_size(.8)
        # Write names of zodiacal constellations between circles 4 and 5
        for i, item in enumerate(text[language]["zodiacal_constellations"]):
            i += 1
            name = "{} {}".format(item['name'], item['symbol'])
            print(name)
            context.circular_text(text=name,
                                  centre_x=0, centre_y=0, radius=rt_2,
                                  azimuth=(-15 + 30 * i),
                                  spacing=1, size=.8)

        # Between circles 5 and 10, display calendars for 2025... 2029 ?
        sp_d = r_5

        # UTC !
        # 2024-03-20	03:06:21	2024-06-20	20:50:56	2024-09-22	12:43:36	2024-12-21	09:20:30
        # 2025-03-20	09:01:25	2025-06-21	02:42:11	2025-09-22	18:19:16	2025-12-21	15:03:01
        # 2026-03-20	14:45:57	2026-06-21	08:24:30	2026-09-23	00:05:13	2026-12-21	20:50:14
        # 2027-03-20	20:24:41	2027-06-21	14:10:50	2027-09-23	06:01:43	2027-12-22	02:42:10
        # 2028-03-20	02:17:08	2028-06-20	20:02:00	2028-09-22	11:45:18	2028-12-21	08:19:40
        # 2029-03-20	08:01:59	2029-06-21	01:48:18	2029-09-22	17:38:30	2029-12-21	14:14:06
        # 2030-03-20	13:52:06	2030-06-21	07:31:19	2030-09-22	23:26:53	2030-12-21	20:09:38


        # UTC! NEED TO ADJUST TO TZ
        equi_solt = (
            #"2024-03-20 03:06:21", "2024-06-20 20:50:56",
                                                          "2024-09-22 12:43:36", "2024-12-21 09:20:30",
            "2025-03-20 09:01:25", "2025-06-21 02:42:11", "2025-09-22 18:19:16", "2025-12-21 15:03:01",
            "2026-03-20 14:45:57", "2026-06-21 08:24:30", "2026-09-23 00:05:13", "2026-12-21 20:50:14",
            "2027-03-20 20:24:41", "2027-06-21 14:10:50", "2027-09-23 06:01:43", "2027-12-22 02:42:10",
            "2028-03-20 02:17:08", "2028-06-20 20:02:00", "2028-09-22 11:45:18", "2028-12-21 08:19:40",
            "2029-03-20 08:01:59", "2029-06-21 01:48:18", "2029-09-22 17:38:30", "2029-12-21 14:14:06",
            "2030-03-20 13:52:06", "2030-06-21 07:31:19", "2030-09-22 23:26:53", "2030-12-21 20:09:38",
        )

        sp_ext = r_5

        last_month = [0] * 13   # idx 0 not used
        for trim in range(-1, 17):

            sp_int12 = sp_ext - dr / 2
            sp_int23 = sp_ext - dr * 2 / 3
            sp_int = sp_ext - dr

            q1 = datetime.fromisoformat(equi_solt[trim + 1]) + timedelta(hours=settings['timezone'])
            q2 = datetime.fromisoformat(equi_solt[trim + 2]) + timedelta(hours=settings['timezone'])
            ts1 = q1.timestamp()
            ts2 = q2.timestamp()

            t = q1.timetuple()
            assert( t[3:] != (0, 0, 0) )
            day_inc = datetime(t[0], t[1], t[2] + 1, 0, 0, 0)
            tsf = day_inc.timestamp()

            inc = pi/2 / ((ts2 - ts1) / 24 / 3600)
            ang =  - pi /2 + (pi/2) * trim + inc * (tsf - ts1)  / 24 / 3600

            if (trim // 2) % 2 == 0 or trim < 0:
                cy = 0
            else:
                cy = - dr/2

            #print(inc, sp_ext, cy)

            while True:
                day_inc += timedelta(days=1)
                t = day_inc.timetuple()
                if t[2] == 1:
                    line_width = 2
                    sp = sp_int
                    print('ang month', t[1], ang)
                    last_month[t[1]] = ang % (2*pi)
                elif t[2] % 10 == 1:
                    line_width = 1
                    sp = sp_int
                elif t[2] % 5 == 1:
                    line_width = 1
                    sp = sp_int23
                else:
                    line_width = 1
                    sp = sp_int12

                if (0 <= trim or t[1] == 12):
                    context.begin_path()
                    context.move_to(x=sp_ext * cos(ang), y=cy-sp_ext * sin(ang))
                    context.line_to(x=sp * cos(ang), y=cy-sp * sin(ang))
                    context.stroke(line_width = line_width)
                if q2 < day_inc or (trim == 16 and t[1] == 1 and t[2] == 1):
                    break
                ang += inc

            if 0 < trim and trim % 2 == 1:
                sp_ext -= dr/2

        for x in range(5):
            sp_d -= dr/2
            context.begin_path()
            context.arc(centre_x=0, centre_y=0 - dr/2, radius=sp_d, arc_from= pi/2, arc_to= 3 * pi/2)
            context.stroke(1)
            if x != 4:
                sp_d -= dr/2
                context.begin_path()
                context.arc(centre_x=0, centre_y=0, radius=sp_d, arc_from= - pi/2, arc_to= pi/2)
                context.stroke(1)

        sp_d -= dr/2
        context.begin_path()
        context.arc(centre_x=0, centre_y=0, radius=sp_d, arc_from = 2*pi - last_month[1] - 0.0027, arc_to = pi/2)
        context.stroke(1)

        context.circular_text(text=str(t[0]), centre_x=0, centre_y=0, radius=sp_ext-0.0004,
                              azimuth= last_month[1] / unit_deg + 2, spacing=1, size=0.5)

        context.circular_text(text=equi_solt[0][:4], centre_x=0, centre_y=0, radius=r_5-0.0004,
                              azimuth= last_month[12] / unit_deg - 2, spacing=1, size=0.5)


        cy = 0
        #sp_ext = sp_d
        sp_int = r_12
        for m in range(12):
            ang = last_month[m+1]
            if 6 <= m :
                print('z', m, ang)
                cy = - dr/2 + 0.0001
                c =  (ang - pi/2)/pi *dr
                sp_int = r_12 - dr/ 2 + c
            context.begin_path()
            context.move_to(x=sp_ext * cos(ang), y=cy -sp_ext * sin(ang))
            context.line_to(x=sp_int * cos(ang), y=cy -sp_int * sin(ang))
            context.stroke(1)

        r_10 -= 0.0002
        # Label names of months
        for mn, (mlen, name) in enumerate(text[language]['months']):
            if mn < 6:
                r = r_10 + dr/2
            else:
                r = r_10 +  (12 - mn) * dr/12

            context.circular_text(text=name, centre_x=0, centre_y=0, radius=r,
                                  azimuth= (last_month[mn+1] + pi / 12)/ unit_deg, spacing=1, size=0.8)





        # Between circles 5 and 10, display calendars for 1394 and 1974

        # The Tuckerman tables provide the longitude of the Sun along the ecliptic on any given day of the year.
        # We produce functions which interpolate the tabulated longitudes, so that we can look up the longitude
        # of the Sun at any moment in time.

        #x_1974 = []
        #y_1974 = []
        #with open("raw_data/tuckerman.dat", "rt") as f_in:
        #    for line in f_in:
        #        line = line.strip()
        #
        #        # Ignore blank lines and comment lines
        #        if (len(line) == 0) or (line[0] == '#'):
        #            continue
        #
        #        # Split line into words
        #        columns = [float(i) for i in line.split()]
        #
        #        x_1974.append(calendar.julian_day(year=1974, month=int(columns[0]), day=int(columns[1]),
        #                                          hour=12, minute=0, sec=0))
        #
        #        y_1974.append(30 * unit_deg * (columns[6] - 1) + columns[7] * unit_deg)

        # Use scipy to do linear interpolation between the supplied data
        #theta_1974 = scipy.interpolate.interp1d(x=x_1974, y=y_1974, kind='linear')

        # Mark 365 days around calendar using the solar longitude data we have.
        # Write numbers on the 10th, 20th and last day of each month

        #rt_1 = (r_6 + r_7) / 2  # Radius of text for the 1394 calendar
        #rt_2 = (r_8 + r_9) / 2  # Radius of text for the 1974 calendar

        #prev_theta = 30 * unit_deg * (10 - 1) + 9.4 * unit_deg
        #
        #with open("raw_data/tuckerman.dat") as f_in:
        #    for line in f_in:
        #        line = line.strip()
        #
        #        # Ignore blank lines and comment lines
        #        if (len(line) == 0) or (line[0] == "#"):
        #            continue

                #m, d, interval, last, z1394, a1394, z1974, a1974 = [float(i) for i in line.split()]

                # *** Calendar for 1974 ***

                # Work out azimuth of given date in 1974 calendar
                #theta = 30 * unit_deg * (z1974 - 1) + a1974 * unit_deg

                # Interpolate interval into number of days since last data point in table (normally five)
                #if prev_theta > theta:
                #    prev_theta = prev_theta - unit_rev
                #for i in arange(0, interval - 0.1):
                #    theta_day = prev_theta + (theta - prev_theta) * (i + 1) / interval
                #    context.begin_path()
                #    context.move_to(x=r_5 * cos(theta_day), y=-r_5 * sin(theta_day))  # r_7 -> r_5
                #    context.line_to(x=r_8 * cos(theta_day), y=-r_8 * sin(theta_day))
                #    context.stroke()
                #prev_theta = theta

                # Draw a marker line on calendar. Month ends get longer markers
                #if last:
                #    context.begin_path()
                #    context.move_to(x=r_8 * cos(theta), y=-r_8 * sin(theta))
                #    context.line_to(x=r_10 * cos(theta), y=-r_10 * sin(theta))
                #    context.stroke()
                #else:
                #    context.begin_path()
                #    context.move_to(x=r_8 * cos(theta), y=-r_8 * sin(theta))
                #    context.line_to(x=r_9 * cos(theta), y=-r_9 * sin(theta))
                #    context.stroke()

                # Label 10th and 20th day of month, and last day of month
                #if ((d % 10) == 0) or (d > 26):
                #    context.set_font_size(1.0)
                #    theta2 = theta - 0.2 * unit_deg
                #    context.text(text="{:.0f}".format(floor(d / 10)),
                #                 x=rt_2 * cos(theta2), y=-rt_2 * sin(theta2),
                #                 h_align=1, v_align=0, gap=0, rotation=-theta - 90 * unit_deg)
                #    theta2 = theta + 0.2 * unit_deg
                #    context.text(text="{:.0f}".format(d % 10),
                #                 x=rt_2 * cos(theta2), y=-rt_2 * sin(theta2),
                #                 h_align=-1, v_align=0, gap=0, rotation=-theta - 90 * unit_deg)

                # *** Calendar for 1394 ***

                # Work out azimuth of given date in 1394 calendar
                #if settings['astrolabe_type'] == 'full':
                #    theta = 30 * unit_deg * (z1394 - 1) + a1394 * unit_deg
                #    if last:
                #        context.begin_path()
                #        context.move_to(x=r_5 * cos(theta), y=-r_5 * sin(theta))
                #        context.line_to(x=r_7 * cos(theta), y=-r_7 * sin(theta))
                #        context.stroke()
                #    else:
                #        context.begin_path()
                #        context.move_to(x=r_6 * cos(theta), y=-r_6 * sin(theta))
                #        context.line_to(x=r_7 * cos(theta), y=-r_7 * sin(theta))
                #        context.stroke()
                #
                #    # Label 10th and 20th day of month, and last day of month
                #    if ((d % 10) == 0) or (d > 26):
                #        context.set_font_size(0.75)
                #        theta2 = theta - 0.2 * unit_deg
                #        context.text(text="{:.0f}".format(d / 10),
                #                     x=rt_1 * cos(theta2), y=-rt_1 * sin(theta2),
                #                     h_align=1, v_align=0, gap=0, rotation=-theta - 90 * unit_deg)
                #        theta2 = theta + 0.2 * unit_deg
                #        context.text(text="{:.0f}".format(d % 10),
                #                     x=rt_1 * cos(theta2), y=-rt_1 * sin(theta2),
                #                     h_align=-1, v_align=0, gap=0, rotation=-theta - 90 * unit_deg)


        ## Add dates of saints days between circles 10 and 12
        #context.set_font_size(1.0)
        #with open("raw_data/saints_days.dat") as f_in:
        #    for line in f_in:
        #        line = line.strip()
        #
        #        # Ignore blank lines and comment lines
        #        if (len(line) == 0) or (line[0] == "#"):
        #            continue
        #
        #        d, m, name = line.split()
        #
        #        day_week = floor(calendar.julian_day(year=1974, month=int(m), day=int(d), hour=12, minute=0, sec=0) -
        #                         calendar.julian_day(year=1974, month=1, day=1, hour=12, minute=0, sec=0)) % 7
        #        sunday_letter = "abcdefg"[day_week:day_week + 1]
        #        theta = theta_1974(calendar.julian_day(year=1974, month=int(m), day=int(d), hour=12, minute=0, sec=0))
        #
        #        #context.circular_text(text=name, centre_x=0, centre_y=0, radius=r_10 * 0.65 + r_11 * 0.35,
        #        #                      azimuth=theta / unit_deg,
        #        #                      spacing=1, size=1)
        #        #context.circular_text(text=sunday_letter, centre_x=0, centre_y=0, radius=r_11 * 0.65 + r_12 * 0.35,
        #        #                      azimuth=theta / unit_deg,
        #        #                      spacing=1, size=1)




        # Shadow scale in middle of astrolabe
        if settings['astrolabe_type'] == 'full':
            context.begin_path()

            # Draw horizontal radial line labelled Occidens
            theta_a = 0 * unit_deg
            context.move_to(x=r_12 * cos(theta_a), y=-r_12 * sin(theta_a))
            context.line_to(x=r_13 * cos(theta_a), y=-r_13 * sin(theta_a))

            # Radial line between RECTA and VERSA
            theta_b = - 45 * unit_deg
            context.move_to(x=r_12 * cos(theta_b), y=-r_12 * sin(theta_b))
            context.line_to(x=r_13 * cos(theta_b), y=-r_13 * sin(theta_b))

            # Radial line between UMBRA and UMBRA
            theta_c = -135 * unit_deg
            context.move_to(x=r_12 * cos(theta_c), y=-r_12 * sin(theta_c))
            context.line_to(x=r_13 * cos(theta_c), y=-r_13 * sin(theta_c))

            # Draw horizontal radial line labelled Oriens
            theta_d = 180 * unit_deg
            context.move_to(x=r_12 * cos(theta_d), y=-r_12 * sin(theta_d))
            context.line_to(x=r_13 * cos(theta_d), y=-r_13 * sin(theta_d))

            # Vertical line at right edge of shadow scale
            context.move_to(x=r_12 * cos(theta_b), y=-r_12 * sin(theta_b))
            context.line_to(x=r_12 * cos(theta_b), y=0)

            # Horizontal line along bottom of shadow scale
            context.move_to(x=r_12 * cos(theta_b), y=-r_12 * sin(theta_b))
            context.line_to(x=r_12 * cos(theta_c), y=-r_12 * sin(theta_c))

            # Vertical line at left edge of shadow scale
            context.move_to(x=r_12 * cos(theta_c), y=-r_12 * sin(theta_c))
            context.line_to(x=r_12 * cos(theta_c), y=0)

            # Central vertical line down middle of shadow scale
            context.move_to(x=0, y=-r_12 * sin(theta_c))
            context.line_to(x=0, y=r_13)
            context.move_to(x=0, y=-r_12)
            context.line_to(x=0, y=-r_13)

            rs1 = r_12 - 0.75 * d_12 / 2  # Radius of corners of fine shadow scale

            rs2 = rs1 - 0.75 * d_12  # Radius of corners of coarse shadow scale

            # Draw horizontal and vertical sides of the fine and coarse shadow scales
            context.move_to(x=rs1 * cos(theta_b), y=-rs1 * sin(theta_b))
            context.line_to(x=rs1 * cos(theta_b), y=0)
            context.move_to(x=rs1 * cos(theta_b), y=-rs1 * sin(theta_b))
            context.line_to(x=rs1 * cos(theta_c), y=-rs1 * sin(theta_c))
            context.move_to(x=rs1 * cos(theta_c), y=-rs1 * sin(theta_c))
            context.line_to(x=rs1 * cos(theta_c), y=0)
            context.move_to(x=rs2 * cos(theta_b), y=-rs2 * sin(theta_b))
            context.line_to(x=rs2 * cos(theta_b), y=0)
            context.move_to(x=rs2 * cos(theta_b), y=-rs2 * sin(theta_b))
            context.line_to(x=rs2 * cos(theta_c), y=-rs2 * sin(theta_c))
            context.move_to(x=rs2 * cos(theta_c), y=-rs2 * sin(theta_c))
            context.line_to(x=rs2 * cos(theta_c), y=0)

            context.stroke()

            # Write the UMBRA and VERSA labels on the shadow scale
            context.set_font_size(0.64)
            context.text(text="UMBRA", x=-1 * unit_mm, y=-rs2 * sin(theta_c), h_align=1, v_align=-1, gap=0.7 * unit_mm,
                         rotation=0)
            context.text(text="UMBRA", x=rs2 * cos(theta_c), y=unit_mm, h_align=-1, v_align=-1, gap=0.7 * unit_mm,
                         rotation=pi / 2)
            context.text(text="RECTA", x=1 * unit_mm, y=-rs2 * sin(theta_c), h_align=-1, v_align=-1, gap=0.7 * unit_mm,
                         rotation=0)
            context.text(text="VERSA", x=rs2 * cos(theta_b), y=unit_mm, h_align=1, v_align=-1, gap=0.7 * unit_mm,
                         rotation=-pi / 2)
            context.text(text="ORIENS", x=-r_12 * 0.95, y=0, h_align=-1, v_align=-1, gap=0.8 * unit_mm, rotation=0)
            context.text(text="OCCIDENS", x=r_12 * 0.95, y=0, h_align=1, v_align=-1, gap=0.8 * unit_mm, rotation=0)

            #r_label = (rs1 + rs2) / 2
            r_label = rs1 - d_12 * 1.5
            offset = 0 # 5 * unit_deg

            # Divisions of scale on shadow scale
            q = 90 * unit_deg
            for i in range(1, 12):
                # Decide how long to make this tick
                rs = rs2 if (i % 4 == 0) else rs1

                # Draw a tick on the shadow scale (right side)
                theta = -atan(i / 12)
                context.begin_path()
                context.move_to(x=rs * cos(theta_b), y=-rs * cos(theta_b) * tan(theta))
                context.line_to(x=r_12 * cos(theta_b), y=-r_12 * cos(theta_b) * tan(theta))
                context.stroke()

                # Label every fourth tick
                if i % 4 == 0:
                    context.text(text="{:d}".format(i),
                                 x=r_label * cos(theta_b), y=-r_label * cos(theta_b) * tan(theta - offset),
                                 h_align=0, v_align=0, gap=0, rotation=-q - theta)

                # Draw a tick on the shadow scale (bottom right)
                theta = -atan(12 / i)
                context.begin_path()
                context.move_to(x=rs * sin(theta_b) / tan(theta), y=-rs * sin(theta_b))
                context.line_to(x=r_12 * sin(theta_b) / tan(theta), y=-r_12 * sin(theta_b))
                context.stroke()

                # Label every fourth tick
                if i % 4 == 0:
                    context.text(text="{:d}".format(i),
                                 x=r_label * sin(theta_b) / tan(theta - offset), y=-r_label * sin(theta_b),
                                 h_align=0, v_align=0, gap=0, rotation=-q - theta)

                # Draw a tick on the shadow scale (bottom left)
                theta = -2 * q - theta
                context.begin_path()
                context.move_to(x=rs * sin(theta_b) / tan(theta), y=-rs * sin(theta_b))
                context.line_to(x=r_12 * sin(theta_b) / tan(theta), y=-r_12 * sin(theta_b))
                context.stroke()

                # Label every fourth tick
                if i % 4 == 0:
                    context.text(text="{:d}".format(i),
                                 x=r_label * sin(theta_b) / tan(theta + offset), y=-r_label * sin(theta_b),
                                 h_align=0, v_align=0, gap=0, rotation=-q - theta)

                # Draw a tick on the shadow scale (left side)
                theta = -2 * q + atan(i / 12)
                context.begin_path()
                context.move_to(x=rs * cos(theta_c), y=-rs * cos(theta_c) * tan(theta))
                context.line_to(x=r_12 * cos(theta_c), y=-r_12 * cos(theta_c) * tan(theta))
                context.stroke()

                # Label every fourth tick
                if i % 4 == 0:
                    context.text(text="{:d}".format(i),
                                 x=r_label * cos(theta_c), y=-r_label * cos(theta_c) * tan(theta + offset),
                                 h_align=0, v_align=0, gap=0, rotation=-q - theta)

            # Add the 12s to the ends of the shadow scale
            #theta = - 45 * unit_deg
            #context.text(text="12", x=r_label * sin(theta_b) / tan(theta - offset), y=-r_label * sin(theta_b),
            #             h_align=0, v_align=0, gap=0, rotation=-pi / 4)
            #
            #theta = -135 * unit_deg
            #context.text(text="12", x=r_label * sin(theta_b) / tan(theta + offset), y=-r_label * sin(theta_b),
            #             h_align=0, v_align=0, gap=0, rotation=pi / 4)

            # Unequal hours scale -- the maths behind this is explained in
            # http://adsabs.harvard.edu/abs/1975JBAA...86...18E

            # First draw innermost circle, which touches centre of astrolabe and the top of the unequal hours scale

            context.set_color(color=theme['blue'])

            context.begin_path()
            context.circle(centre_x=0, centre_y=-r_12 / 2, radius=r_12 / 2)
            context.stroke()

            # Now draw arcs for the hours 1 to 11
            for h in range(1, 12):
                theta = 15 * unit_deg * h
                # Vertical position of the centre of the arc
                y_centre = r_12 * cos(theta) / 2 + r_12 * sin(theta) / 2 * tan(theta)

                # Size of arc
                arc_end = atan2(r_12 * sin(theta), r_12 * cos(theta) / 2 - r_12 * sin(theta) / 2 * tan(theta))

                context.begin_path()
                context.arc(centre_x=0, centre_y=-y_centre, radius=y_centre,
                            arc_from=arc_end - pi / 2, arc_to=-arc_end - pi / 2)
                context.stroke()

                context.text(text="{:d}".format(6 - abs(h - 6)),
                     x=r_12 * cos(theta), y=-r_12 * sin(theta),
                     h_align=0, v_align=2, gap=0, rotation=-theta - pi/2)

        # Finish up
        #context.set_color(color=theme['text'])
        #context.circular_text(text=text[language]['copyright'], centre_x=0, centre_y=0, radius=r_12 - 2 * unit_mm,
        #                      azimuth=270, spacing=1, size=0.7)


        # B(d) = 2 pi ( d - 81) / 365
        # Tc(d) = 7.53 cos(B) + 1.5 sin(B)
        #
        def time_equation(ang):
            return 7.678 * sin(ang + 1.374) - 9.87 * sin(2 * ang)

        r = r_12 * 11.0 / 20     # between r_12 amd 1/4 r_12
        amp = r_12 * 9.0 / 20 / 16    # 16 mn -> 5/8 r_12

        context.set_color(color=theme['red'])
        context.set_line_style(dotted=True)
        context.begin_path()
        context.circle(centre_x=0, centre_y=0, radius=r - amp * 15)
        context.circle(centre_x=0, centre_y=0, radius=r)
        context.circle(centre_x=0, centre_y=0, radius=r + amp * 15)
        context.stroke(line_width=.5)

        context.set_line_style(dotted=False)

        context.begin_path()
        m0 = time_equation(0) * amp + r
        context.move_to(x=m0 * cos(0), y=m0 * sin(0))
        for ang in range(360):
            m = time_equation(ang*unit_deg) * amp + r
            #print(ang, ang*unit_deg, "->", time_equation(ang / unit_deg), m)
            context.line_to(x=m * cos(ang*unit_deg), y=-m * sin(ang*unit_deg))

        context.line_to(x=m0 * cos(0), y=m0 * sin(0))
        context.stroke(line_width=1.5)



# Do it right away if we're run as a script
if __name__ == "__main__":
    # Fetch command line arguments passed to us
    arguments = fetch_command_line_arguments(default_filename=MotherBack().default_filename())

    # Render the back of the mother
    MotherBack(settings={
        'latitude': arguments['latitude'],
        'language': 'en'
    }).render_to_file(
        filename=arguments['filename'],
        img_format=arguments['img_format']
    )
