#!/usr/bin/python3
# climate.py
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
Render the climate of the astrolabe.
"""

from math import pi, sin, tan, cos, atan2, hypot, acos
from typing import Dict

from constants import unit_deg, unit_cm, unit_mm, inclination_ecliptic, centre_scaling, r_1, d_12, tab_size
from graphics_context import BaseComponent, GraphicsContext
from numpy import arange
from settings import fetch_command_line_arguments
from text import text
from themes import themes


class Climate(BaseComponent):
    """
    Render the climate of the astrolabe.
    """

    def default_filename(self) -> str:
        """
        Return the default filename to use when saving this component.
        """
        return "climate"

    def bounding_box(self, settings: dict) -> Dict[str, float]:
        """
        Return the bounding box of the canvas area used by this component.

        :param settings:
            A dictionary of settings required by the renderer.
        :return:
            Dictionary with the elements 'x_min', 'x_max', 'y_min' and 'y_max' set
        """

        r_outer: float = r_1 - d_12 * 2.5

        return {
            'x_min': -r_outer,
            'x_max': r_outer,
            'y_min': -r_outer,
            'y_max': r_outer
        }

    def do_rendering(self, settings: dict, context: GraphicsContext) -> None:
        """
        This method is required to actually render this item.

        :param settings:
            A dictionary of settings required by the renderer.
        :param context:
            A GraphicsContext object to use for drawing
        :return:
            None
        """

        is_southern: bool = settings['latitude'] < 0
        latitude: float = abs(settings['latitude'])
        placename: float = settings['placename']
        language: str = settings['language']
        theme: Dict[str, Tuple[float, float, float, float]] = themes[settings['theme']]

        context.set_color(color=theme['lines'])

        # Define the radii of all the concentric circles drawn on front of mother

        # The radius of the tab at the top of climate, relative to the centre of the astrolabe
        r_tab: float = r_1 - d_12 * 2.5 - unit_mm

        # Outer radius of climate
        r_2: float = r_1 - d_12 * 3 - unit_mm

        # Radius of central hole
        r_3: float = d_12 * centre_scaling

        # Radius of the line denoting the equator
        r_4: float = r_2 * tan((90 - inclination_ecliptic) / 2 * unit_deg)

        # Radius of the line denoting the tropic of Cancer
        r_5: float = r_4 * tan((90 - inclination_ecliptic) / 2 * unit_deg)

        # Draw the outer edge of climate, and the central hole, and use these to create a clipping region
        context.begin_path()
        context.circle(centre_x=0, centre_y=0, radius=r_2)
        context.begin_sub_path()
        context.circle(centre_x=0, centre_y=0, radius=r_3)
        context.stroke()
        #context.rectangle(x0= -0.8*unit_cm, y0= r_2-2.1*unit_cm, x1=.8*unit_cm , y1=r_2 - 1.65*unit_cm)
        #context.rectangle(x0= -0.8*unit_cm, y0= r_2-1.55*unit_cm, x1=.8*unit_cm , y1=r_2 - 1.05*unit_cm)
        context.clip()

        # Make the tab at the top of the climate
        context.begin_path()
        context.arc(centre_x=0, centre_y=0, radius=r_tab,
                    arc_from=-tab_size - pi / 2, arc_to=tab_size - pi / 2)
        context.move_to(x=r_tab * sin(tab_size), y=-r_tab * cos(tab_size))
        context.line_to(x=r_2 * sin(tab_size), y=-r_2 * cos(tab_size))
        context.move_to(x=-r_tab * sin(tab_size), y=-r_tab * cos(tab_size))
        context.line_to(x=-r_2 * sin(tab_size), y=-r_2 * cos(tab_size))
        context.stroke()

        context.set_color(theme['red'])
        # Draw the equator
        context.begin_path()
        context.circle(centre_x=0, centre_y=0, radius=r_4)
        context.stroke()

        # Draw the tropic of Cancer
        context.begin_path()
        context.circle(centre_x=0, centre_y=0, radius=r_5)
        context.stroke()

        context.begin_path()
        context.move_to(x=-r_2, y=0)
        context.line_to(x=r_2, y=0)
        context.move_to(x=0, y=r_2 if settings['astrolabe_type'] == 'full' else r_4)
        context.line_to(x=0, y=-r_2)
        context.stroke(line_width=1, dotted=False)

        # The maths involved in drawing the climate is described in this paper:
        # http://adsabs.harvard.edu/abs/1976JBAA...86..125E

        # Draw lines of constant altitude
        altitude: float
        horizon_centre: float = 0
        horizon_radius: float = 0


        for altitude in [-18, -12, -6] + list(range(0,80,2)) + [70, 80]: #+ list(range(5,70,10)
            theta1: float = (-latitude - (90 - altitude)) * unit_deg
            theta2: float = (-latitude + (90 - altitude)) * unit_deg

            x1: float = r_4 * sin(theta1)
            y1: float = r_4 * cos(theta1)
            x2: float = r_4 * sin(theta2)
            y2: float = r_4 * cos(theta2)

            y_a: float = y1 * (r_4 / (r_4 - x1))
            y_b: float = y2 * (r_4 / (r_4 - x2))

            # Record centre and radius of the arc denoting the horizon
            if altitude == 0:
                horizon_centre = (y_a + y_b) / 2
                horizon_radius = (y_b - y_a) / 2

            context.set_font_style(bold=False)
            context.set_color(theme['text'])
            context.set_font_size(0.60)

            #if y_b < r_2:
            #if False:
            #    if (altitude % 10) == 0:
            #        context.text(text="{:.0f}".format(float(altitude)),
            #                     x=0, y=-y_b, h_align=0, v_align=1, gap=0, rotation=0)
            #else:
            if (0 < altitude) and (altitude % 5 == 0):
                r: float = (y_b - y_a) / 2 * 0.98 - .5 * unit_mm
                y: float = (y_a + y_b) / 2
                print("XXXXX",  (r ** 2 + y ** 2 - r_2 ** 2) ,  y_b - y_a , y_a + y_b, "=",
                    (r ** 2 + y ** 2 - r_2 ** 2) / (2 * ((y_b - y_a) / 2) * ((y_a + y_b) / 2)))
                #start: float = 180 * unit_deg - acos(
                #    (r ** 2 + y ** 2 - r_2 ** 2) / (2 * ((y_b - y_a) / 2) * ((y_a + y_b) / 2)))
                start: float = (108  - altitude*.8 - (20 if altitude == 80 else 0)) * unit_deg
                end: float = -start
                context.text(text="{:.0f}".format(float(altitude)),
                             x=r * sin(start + (r_2 / r) * 3 * unit_deg),
                             y=-(y_a + y_b) / 2 - r * cos(start + (r_2 / r) * 3 * unit_deg),
                             h_align=0, v_align=- altitude / 100,
                             gap=0,
                             rotation=180 * unit_deg + (start + (r_2 / r) * 3 * unit_deg))
                context.text(text="{:.0f}".format(float(altitude)),
                             x=r * sin(end - (r_2 / r) * 2 * unit_deg),
                             y=-(y_a + y_b) / 2 - r * cos(end - (r_2 / r) * 3 * unit_deg),
                             h_align=0, v_align=-0.2,
                             gap=0,
                             rotation=180 * unit_deg + (end - (r_2 / r) * 3 * unit_deg))

            context.begin_path()
            context.circle(centre_x=0, centre_y=-(y_a + y_b) / 2, radius=(y_b - y_a) / 2)
            context.stroke(dotted=(altitude < 0),
                           line_width= .3 + .6 * int(altitude == 0)
                                + .4 * int(altitude % 5 == 0)
                                + .3 * int(altitude % 10 == 0)
                                + .4 * int(altitude % 30 == 0),
                           color=theme['text'] if altitude >= 0 else theme['alt_az'])

            if altitude <= 0:
                context.circular_text(text=text[language]['twilight'][altitude],
                                    centre_x=0, centre_y=-(y_a + y_b) / 2,
                                    radius=(y_b - y_a) / 2 + 1.45 * unit_mm,
                                    azimuth=270 + altitude*1.3  + 1.2*latitude,
                                    spacing=1, size=0.5)

                                    # 50: azimuth=330 + altitude*1.3,  #+ 1.2*latitude
                                    # 40: azimuth=318 + altitude*1.3,  #+ 1.2*latitude

                print("LAT", latitude, 282 - latitude/2 + altitude*0.73)

        # Find coordinates of P
        theta: float = -latitude * unit_deg
        p_x: float = r_4 * sin(theta)
        p_y: float = r_4 * cos(theta)

        # Find coordinates of Z
        z_x: float = 0
        z_y: float = p_y / (r_4 - p_x) * r_4

        # Find midpoint between Z and H
        zh_x: float = -r_4 / 2
        zh_y: float = z_y / 2

        # Find bearing of T from ZH (clockwise from right-going axis)
        theta: float = atan2(z_x - (-r_4), z_y)

        # Find coordinates of T
        t_x: float = 0
        t_y: float = zh_y + zh_x * tan(theta)


        print("P:{:.2f} {:.2f}    Z:{:.2f} {:.2f}    T:{:.2f} {:.2f}".format(p_x/unit_cm, p_y/unit_cm, z_x/unit_cm, z_y/unit_cm, t_x/unit_cm, t_y/unit_cm))

        # Draw lines of constant azimuth. We draw 16 arcs at 11.25 degree intervals, which cut through the zenith
        # and meet the horizon in two opposite compass bearings. For this reason we only draw half as many arcs as
        # there are compass bearings

        step_size: float = 11.25 * unit_deg
        for azimuth_step in range(1, 16):
            azimuth: float = -90 * unit_deg + step_size * azimuth_step

            t_x: float = (z_y - t_y) * tan(azimuth)

            # Radius of arc of constant azimuth
            t_r: float = hypot(t_x, t_y - z_y)
            t_r2: float = t_r * 1.0

            t_hc: float = hypot(t_x, t_y - horizon_centre)  # Distance from T to centre of horizon
            theta: float = acos((t_r ** 2 + t_hc ** 2 - horizon_radius ** 2) / (2 * t_r * t_hc))
            phi: float = atan2(t_x, horizon_centre - t_y)
            start: float = -phi - theta
            end: float = -phi + theta

            t_c: float = hypot(t_x, t_y)  # Distance from T to centre of the astrolabe
            arg: float = (t_r ** 2 + t_c ** 2 - r_2 ** 2) / (2 * t_r * t_c)
            if (arg >= 1) or (arg <= -1):
                start2: float = start
                end2: float = end
            else:
                theta: float = acos((t_r ** 2 + t_c ** 2 - r_2 ** 2) / (2 * t_r * t_c))
                phi: float = atan2(t_x, -t_y)
                start2: float = -phi - theta
                end2: float = -phi + theta

            context.begin_path()
            context.arc(centre_x=t_x, centre_y=-t_y, radius=t_r,
                        arc_from=max(start, start2) - pi / 2, arc_to=min(end, end2) - pi / 2)
            context.stroke(line_width=0.5 if azimuth else 1,
                           color=theme['alt_az'])

            # Compass direction for the start and end of the line of constant azimuth. Each line of constant azimuth
            # meets the horizon at two opposite points, with opposite compass directions.
            if azimuth_step % 2 == 0:

                def float2text(f):
                    if f % 1 == 0:
                        return str(int(f)) + '°'
                    else:
                        return str(f) + '°'

                direction_start: str = text[language]['directions'][azimuth_step // 2]
                direction_end: str = text[language]['directions'][azimuth_step // 2 + 8]

                # In southern hemisphere, invert directions
                if is_southern:
                    direction_start, direction_end = (direction_end, direction_start)

                context.set_color(theme['blue'])
                if hypot(t_x + t_r * sin(end), t_y + t_r * cos(end)) < 0.9 * r_2:
                    context.set_font_size(0.90)
                    context.set_font_style(bold=True)
                    context.text(text=direction_start,
                                 x=t_x + t_r * sin(end), y=-t_y - t_r * cos(end),
                                 h_align=0, v_align=1, gap=unit_mm,
                                 rotation=end - 90 * unit_deg)

                    context.set_font_size(0.50)
                    context.set_font_style(bold=False)
                    context.text(text=float2text(270 - azimuth/unit_deg),
                                 x=t_x + t_r * sin(end), y=-t_y - t_r * cos(end),
                                 h_align=0, v_align=-1, gap=unit_mm,
                                 rotation=end - 90 * unit_deg)
                else:
                    context.set_font_size(0.90)
                    context.set_font_style(bold=True)
                    context.text(text=direction_start,
                                 x=t_x + t_r * sin(min(end, end2) - (r_2 / t_r) * 2.1 * unit_deg),
                                 y=-t_y - t_r * cos(min(end, end2) - (r_2 / t_r) * 2.1 * unit_deg),
                                 h_align=0, v_align=0, gap=0,
                                 rotation=90 * unit_deg + min(end, end2) - (r_2 / t_r) * 8 * unit_deg)
                    context.set_font_size(0.50)
                    context.set_font_style(bold=False)
                    context.text(text=float2text(270 - azimuth/unit_deg),
                                 x=t_x + t_r * sin(min(end, end2) - (r_2 / t_r) * 4 * unit_deg),
                                 y=-t_y - t_r * cos(min(end, end2) - (r_2 / t_r) * 4 * unit_deg),
                                 h_align=0, v_align=0, gap=0,
                                 rotation=90 * unit_deg + min(end, end2) - (r_2 / t_r) * 8 * unit_deg)

                if hypot(t_x + t_r * sin(start), t_y + t_r * cos(start)) < 0.9 * r_2:
                    context.set_font_size(0.90)
                    context.set_font_style(bold=True)
                    context.text(text=direction_end,
                                 x=t_x + t_r * sin(start),
                                 y=-t_y - t_r * cos(start),
                                 h_align=0, v_align=1, gap=unit_mm,
                                 rotation=90 * unit_deg + start)

                    context.set_font_size(0.50)
                    context.set_font_style(bold=False)
                    context.text(text=float2text(90 - azimuth/unit_deg),
                                 x=t_x + t_r * sin(start),
                                 y=-t_y - t_r * cos(start),
                                 h_align=0, v_align=-1, gap=unit_mm,
                                 rotation=90 * unit_deg + start)
                else:
                    context.set_font_size(0.90)
                    context.set_font_style(bold=True)
                    context.text(text=direction_end,
                                 x=t_x + t_r * sin(max(start, start2) + (r_2 / t_r) * 2.1 * unit_deg),
                                 y=-t_y - t_r * cos(max(start, start2) + (r_2 / t_r) * 2.1 * unit_deg),
                                 h_align=0, v_align=0, gap=0,
                                 rotation= -90 * unit_deg + max(start, start2) + (r_2 / t_r) * 8 * unit_deg)
                    context.set_font_size(0.50)
                    context.set_font_style(bold=False)
                    context.text(text=float2text(90 - azimuth/unit_deg),
                                 x=t_x + t_r * sin(max(start, start2) + (r_2 / t_r) * 4 * unit_deg),
                                 y=-t_y - t_r * cos(max(start, start2) + (r_2 / t_r) * 4 * unit_deg),
                                 h_align=0, v_align=0, gap=0,
                                 rotation= -90 * unit_deg + max(start, start2) + (r_2 / t_r) * 8 * unit_deg)

        context.set_color(theme['blue'])
        context.set_font_style(bold=True)
        context.set_font_size(0.90)
        context.text(text="N" if not is_southern else "S",
                     x=0, y=-horizon_centre + horizon_radius,
                     h_align=0, v_align=1, gap=unit_mm, rotation=0)
        context.text(text="S" if not is_southern else "N",
                     x=0, y=-r_2,
                     h_align=0, v_align=1, gap=unit_mm, rotation=0)
        context.set_font_size(0.50)
        context.set_font_style(bold=False)
        context.text(text="0°" if not is_southern else "180°",
                     x=0, y=-horizon_centre + horizon_radius,
                     h_align=0, v_align=-1, gap=unit_mm, rotation=0)
        context.text(text="180°" if not is_southern else "0°",
                     x=0, y=-r_2,
                     h_align=0, v_align=2.7, gap=unit_mm, rotation=0)


        # Subroutine for calculating the azimuthal angle of the lines of the unequal hours
        if settings['astrolabe_type'] == 'full':
            # Subroutine for calculating the azimuthal angle of the lines of the unequal hours
            def theta_unequal_hours(r: float) -> float:
                arg: float = (r ** 2 + horizon_centre ** 2 - horizon_radius ** 2) / (2 * r * horizon_centre)
                if arg <= -1:
                    return 180 * unit_deg
                if arg >= 1:
                    return 0 * unit_deg
                return acos(arg)

            context.set_color(theme['red'])
            # Draw lines of unequal hours in turn
            for h in arange(1, 12, 0.5):
                for r in arange(max(r_5, horizon_radius - horizon_centre), r_2 + 0.05 * unit_mm, 0.5 * unit_mm):
                    r0: float = r
                    r1: float = min(r + 0.5 * unit_mm, r_2)
                    theta0: float = theta_unequal_hours(r0)
                    theta1: float = theta_unequal_hours(r1)
                    psi0: float = theta0 + (360 * unit_deg - 2 * theta0) / 12 * h
                    psi1: float = theta1 + (360 * unit_deg - 2 * theta1) / 12 * h
                    context.begin_path()
                    context.move_to(x=r0 * sin(psi0), y=-r0 * cos(psi0))
                    context.line_to(x=r1 * sin(psi1), y=-r1 * cos(psi1))
                    context.stroke(line_width= .2 if h%1 else 1, dotted=False)

            # Label the unequal hours
            context.set_font_size(1.6)
            r: float = r_2 - 4 * unit_mm
            theta0: float = theta_unequal_hours(r)
            context.set_font_style(bold=False)
            for pos, hr in enumerate(["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII"]):
                psi0: float = theta0 + (360 * unit_deg - 2 * theta0) / 12 * (pos + 0.5)
                psi0 = (psi0 - 180 * unit_deg) * 0.95 + 180 * unit_deg
                context.text(text=hr,
                             x=r * sin(psi0), y=-r * cos(psi0),
                             h_align=0, v_align=0, gap=unit_mm,
                             rotation=180 * unit_deg + psi0)

        # A space to write the owner's name
        if settings['astrolabe_type'] != 'full':
            arc_size: float = 40 * unit_deg
            context.begin_path()
            context.move_to(x=r_2 * sin(arc_size), y=r_2 * cos(arc_size))
            context.arc(centre_x=0, centre_y=0,
                        radius=r_2 - 0.8 * unit_cm,
                        arc_from=90 * unit_deg - arc_size,
                        arc_to=90 * unit_deg + arc_size
                        )
            context.line_to(x=-r_2 * sin(arc_size), y=r_2 * cos(arc_size))
            context.stroke(line_width=1, dotted=False)

            context.circular_text(text="{}:".format(text[language]['name']),
                                  centre_x=0, centre_y=0,
                                  radius=r_2 - 0.4 * unit_cm,
                                  azimuth=238,
                                  spacing=1, size=1.2)

        # Draw horizontal and vertical lines through the middle of the climate
        #context.begin_path()
        #context.rectangle(x0= -.5*unit_cm, y0= r_2-2.3*unit_cm, x1=.5*unit_cm , y1=r_2 - 1.6*unit_cm)
        #context.clip()

        context.begin_path()
        context.move_to(x=-r_2, y=0)
        context.line_to(x=r_2, y=0)
        context.move_to(x=0, y=r_2 if settings['astrolabe_type'] == 'full' else r_4)
        context.line_to(x=0, y=-r_2)
        context.stroke(line_width=1, dotted=False)

        #context.clip_reset()

        # Finish up
        context.set_font_style(bold=False)
        #context.circular_text(text=text[language]['url'],
        #                      centre_x=0, centre_y=0, radius=r_2 - 1.6 * unit_cm,
        #                      azimuth=270, spacing=1, size=0.7)
        #context.circular_text(text=text[language]['copyright'],
        #                      centre_x=0, centre_y=0, radius=r_2 - 1.3 * unit_cm,
        #                      azimuth=270, spacing=1, size=0.7)
        #context.circular_text(text=text[language]['climate_latitude'].format(latitude, "N" if not is_southern else "S"),
        #                      centre_x=0, centre_y=0, radius=r_2 - 1.0 * unit_cm,
        #                      azimuth=270, spacing=1, size=0.7)

        #context.begin_path()
        #context.rectangle(x0= -r_2, y0= -r_2, x1= r_2, y1= r_2)
        #context.begin_sub_path()
        #context.rectangle(x0= -.5*unit_cm, y0= r_2-2.3*unit_cm, x1=.5*unit_cm , y1=r_2 - 1.6*unit_cm)
        ###context.rectangle(x0=-100*unit_cm, y0= -100*unit_cm, x1=-99*unit_cm , y1=99*unit_cm)
        #context.stroke(line_width=1, dotted=False)
        #context.clip()

        context.set_color(color=theme['lines'])
        context.text(text=placename,
                     x= 0 ,
                     y= r_2 - 1.6 * unit_cm,
                     h_align=0, v_align=-0.95, gap=unit_mm,
                     rotation=0)

        context.text(text=text[language]['climate_latitude'].format(latitude, "N" if not is_southern else "S"),
                     x= 0 ,
                     y= r_2 - 1.6 * unit_cm,
                     h_align=0, v_align=0.95, gap=unit_mm,
                     rotation=0)


# Do it right away if we're run as a script
if __name__ == "__main__":
    # Fetch command line arguments passed to us
    arguments = fetch_command_line_arguments(default_filename=Climate().default_filename())

    # Render the climate
    Climate(settings={
        'latitude': arguments['latitude'],
        'language': 'en'
    }).render_to_file(
        filename=arguments['filename'],
        img_format=arguments['img_format']
    )
