# OpenGL-Drawing-Filling

## How To Use The Program
 At the top of the GUI you will be given an option to choose a line drawing algorithm
between DDA and Bresenham. Once you select either of the options the polygon(s) will be
drawn. There is a button next to the options to fill the polygons, you can continuously click the
button to fill and unfill the polygons.

You will be given the option to input x-min, x-max, y-min, and y-max values. It is
important to note that a x value of 0 is located at the very left of the window. These values
dictate the region in which the polygons will be visible

You will have the option to choose which polygon you want to manipulate ranging from 0
to n - 1. With the first polygon in the input file having an ID of 0 and the last an ID of n - 1.

You can input an X and Y value which will translate the chosen polygon that many units.
You must click the button below where you input the X and Y values in order to apply the
translation.

You can input any number in which you want to scale the polygon, you must click the
button below where you input the value to apply the scale

There is a slider labled “rotate” moving the slider to left rotates the polygon
counter-clockwise, moving to the right rotates the polygon clockwise. You must click the button
labeled “rotate” to apply the rotation to the polygon.

## Algorithims Implemented

- Scan Fill Algorithm
- Sutherland Hodgeman
- Cohen Sutherland
- Bresenham Line Algorithim
- DDA Line Algorithim
