function [Heading xh yh] = clcHeading(Hxyz, Pitch, Roll, headingGT180, plotMagnetics)
% Function: Calcualates heading from magnetometer data.  Applies correction
% based on Roll/Pitch data.
%
% Parameters:
% Roll: array of bank angles in radians(!)
% Pitch: array of elevation angles in radians
% Hxyz(:,1): array of magnetometer x-axis data in gauss
% Hxyz(:,2): array of magnetometer y-axis data in gauss
% Hxyz(:,3): array of magnetometer z-axis data in gauss
% headingGT180: true if heading should be in the range 0 to +/-180, else
%               false if in range 0-360.
%               Default: false
% plotMagnetics: true if magnetics should be plotted, else false
%               Default: false
% Returns:
% A three element structure containing the elements heading, corrected magnetic
% x and corrected magnetic y values.  Each element is an array with a size 
% equal to the number of rows in the original data file.
%
% Initial concept based on heading paper from Honeywell.  See:
% http://www.honeywell.com/sites/servlet/com.merx.npoint.servlets.DocumentServlet?docid=D84A3A7BE-2A47-A0C3-F3BF-E2F7768BD449
% Essentially, taking out the Roll/Pitch rotations, then calculating x and y components.
    Heading= sin(-Roll);  %temporary
    yh=      sin(-Pitch); %temporary 
    xh= Hxyz(:,1).*cos(-Pitch) + ...
       (Hxyz(:,2).*Heading    - Hxyz(:,3).*cos(-Roll)).*yh; % Gauss                
    yh= Hxyz(:,2).*cos(-Roll) + Hxyz(:,3).*Heading;         % Gauss
       
    % NOTE: will always result in +/-180 unless modified by next block.
    % Also, heading is negated so the results adhere to right-hand-rule 
    % conventions.
    Heading = -atan2(yh, xh); 
    
    % determine default parameter values.
    if nargin < 7, plotMagnetics = false; end;
    if nargin < 6, headingGT180 = false; end;
    % Determine if heading results should be in the form 0-360.  
    % This block can be effective when 'wraparound' at +/-180 is present.
    if headingGT180
        for i = 1:length(Heading)
            if Heading(i) < 0 
                Heading(i) = 2 * pi + Heading(i);
            end;
        end;
    end;
        
    % Map out raw and uncorrected magnetic components
    if plotMagnetics
        % plot raw mag values
        figure;
        subplot(2, 1, 1);
        plot(Hxyz(:,1));
        hold;
        plot(Hxyz(:,2), 'g');
        plot(Hxyz(:,3), 'r');
        ylabel('Gauss');
        title('Magnetic Components');
        legend('Mx', 'My', 'Mz');
        
        % plot out heading based on uncorrected mag data
        subplot(2, 1, 2);
        heading_raw = -atan2(Hxyz(:,2), Hxyz(:,1));
        
        if headingGT180
            for i = 1:length(heading_raw)
                if heading_raw(i) < 0 
                    heading_raw(i) = 2 * pi + heading_raw(i);
                end;
            end;
        end;
        
        heading_raw = rad2deg(heading_raw);
        plot(heading_raw);
        title('Uncorrected Heading');
        ylabel('Heading(deg)');
    end;
end