function check = isPAACentroidValid(roomCoordinates, centroidCoordinate)
%ISPAACENTROIDVALID returns check if the centroid inside the defined area (box).
%
%   ISPAACENTROIDVALID(roomCoordinates, centroidCoordinate)
%   Generate a warning is the centroid coordinates centroidCoordinate are
%   outside the area defined in roomCoordinates
%
%   check = ISPAACENTROIDVALID(roomCoordinates, centroidCoordinate)
%   Return false if the centroid coordinates centroidCoordinate are
%   outside the area defined in roomCoordinates
%   
%--------------------------Software Disclaimer-----------------------------
%
% NIST-developed software is provided by NIST as a public service. You may 
% use, copy and distribute copies of the software in any medium, provided 
% that you keep intact this entire notice. You may improve, modify and  
% create derivative works of the software or any portion of the software, 
% and you  may copy and distribute such modifications or works. Modified 
% works should carry a notice stating that you changed the software and  
% should note the date and nature of any such change. Please explicitly  
% acknowledge the National Institute of Standards and Technology as the 
% source of the software.
% 
% NIST-developed software is expressly provided "AS IS." NIST MAKES NO
% WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT OR ARISING BY OPERATION  
% OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT AND 
% DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF 
% THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS 
% WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS  
% REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT 
% NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF 
% THE SOFTWARE.
%
% You are solely responsible for determining the appropriateness of using
% and distributing the software and you assume all risks associated with  
% its use, including but not limited to the risks and costs of program 
% errors, compliance with applicable laws, damage to or loss of data, 
% programs or equipment, and the unavailability or interruption of 
% operation. This software is not intended to be used in any situation  
% where a failure could cause risk of injury or damage to property. The 
% software developed by NIST employees is not subject to copyright 
% protection within the United States.

% Copyright 2019-2020 NIST/CTL (steve.blandino@nist.gov)

%#codegen

roomCoordinates = reshape(permute(reshape(roomCoordinates,size(roomCoordinates,1),3,[]),[1,3,2]),[],3);
minCoord = min(roomCoordinates);
maxCoord = max(roomCoordinates);


if ~all(reshape(centroidCoordinate>= minCoord & centroidCoordinate<= maxCoord, [],1))
    warning OFF BACKTRACE
    warning('PAA outside defined area')
    warning ON BACKTRACE
    check = false;
else
    check = true;
end

end