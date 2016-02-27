function getTriad(h,varargin)
%       'linestyle' - [ {-} | -- | : | -. | none ]
%       'linewidth' - Used to define thickness lines representing axes 
%                     Default value is 0.50
%                     Value must be finite and greater than zero
%       'matrix'    - 4x4 homogenious transform
%       'parent'    - Axes or other hgtransform handle
%       'scale'     - Used to define the length of each axes line
%                     Default value is 1.00
%                     Value(s) must be finite and greater than zero
%                     A scalar value scales the x, y, and z-axis equally
%                     A 3-element array (e.g. [1,2,3]) scales each of the
%                       axes seperately.
%       'tag'       - String describing object
%       'visible'   - [ {on} | off ]
%                     This property hides the visualization, not the
%                       hgtransform object.
%TODO - create get function for triads

warning('This function is still a work in progress.');