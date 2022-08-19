function [Pitch Roll]= clcPitchRoll(Gxyz)
% Function: Calcualates bank angle and elevation from accelerometer data.
%
% Parameters:
% ax: array of accelerometer x-axis data
% ay: array of accelerometer y-axis data
% az: array of accelerometer z-axis data
%
% Returns:
% A two element structure containing the elements bank, elevation 
% specified in radians.  Each element is an array with a size equal 
% to the number of rows in the original accel data array.
%   
% Note: tilt method taken from Kionix.  See following paper for
% details: http://kionix.com/sensors/application-notes.html,
% document AN005.
% Note that the elevation is negated so that 'nose-up' is recorded as a
% positive elevation.  You can see this is required if you look at the
% sign of the accelerometer values when rotation occurs; X goes
% negative immediately upon a nose-up attitude, and Z goes positive.
% atan2 looks at signs to determine quadrant, and negative X/
% positive Z is interpreted as quadrant 4 (nose down).

sGxyz       =  Gxyz.^2;
sAxAx       =  sqrt(sum(sGxyz(:,[2,3]), 2));
Pitch       = -atan2(Gxyz(:,1), sAxAx); % -90 <= elevation <= 90
sAxAx       =  sqrt(sum(sGxyz(:,[1,3]), 2));
Roll        =  atan2(Gxyz(:,2), sAxAx); % -90 <= bank      <= 90