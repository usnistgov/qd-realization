function reflectionLossdB = getReflectionLossApproach2(MaterialLibrary,...
                materialIndex, multipath)
% getReflectionLossApproach2 returns reflection loss in dB. 
% This function is based on the MATLAB WLAN toolbox.

% NIST-developed software is provided by NIST as a public service. You may
% use, copy and distribute copies of the software in any medium, provided
% that you keep intact this entire notice. You may improve,modify and
% create derivative works of the software or any portion of the software,
% and you may copy and distribute such modifications or works. Modified
% works should carry a notice stating that you changed the software and
% should note the date and nature of any such change. Please explicitly
% acknowledge the National Institute of Standards and Technology as the
% source of the software. NIST-developed software is expressly provided
% "AS IS." NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT OR
% ARISING BY OPERATION OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
% WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
% NON-INFRINGEMENT AND DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS
% THAT THE OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE,
% OR THAT ANY DEFECTS WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY
% REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF,
% INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY,
% OR USEFULNESS OF THE SOFTWARE.
%
% You are solely responsible for determining the appropriateness of using
% and distributing the software and you assume all risks associated with
% its use,including but not limited to the risks and costs of program
% errors, compliance with applicable laws, damage to or loss of data,
% programs or equipment, and the unavailability or interruption of
% operation. This software is not intended to be used in any situation
% where a failure could cause risk of injury or damage to property.
% The software developed by NIST employees is not subject to copyright
% protection within the United States.
%
% 2020-2021 NIST/CTL (neeraj.varshney@nist.gov)

% Get reflectances for V and H
reflectionCoefficient = getReflectance(MaterialLibrary,materialIndex,multipath);
% Calculate reflection loss. It is the sum (in dB) of the mean
% of the V and H reflectances for each reflection.
reflectionLoss = prod(sum(reflectionCoefficient)/2);
reflectionLossdB = -20*log10(abs((reflectionLoss))); 
end

function reflectionCoefficient = getReflectance(MaterialLibrary,materialIndex,multipath)
% Calculate angles of incident 
incidentAngle = angleOfIncident(multipath);
orderReflection = multipath(1,1);
reflectionCoefficient = ones(2, orderReflection);
for reflectionOrderIndex = 1:orderReflection
    if reflectionOrderIndex == 1
        reflectionMaterialIndex = materialIndex(1,1);
    elseif  reflectionOrderIndex == 2
        reflectionMaterialIndex = materialIndex(1,2);
    else
        error(strcat('Incident angles are obtained till second order ',...
        'reflection. Thus, order of reflection cannot be considered ',... 
        'higher than second order reflection.'));
    end
    % Use Fresnel equation to derive power reflectivity                
    relativePermittivity = MaterialLibrary.RelativePermittivity(reflectionMaterialIndex); 
    aor = incidentAngle(reflectionOrderIndex);
    B_h =  relativePermittivity - sind(aor)^2;                
    B_v = (relativePermittivity - sind(aor)^2)/relativePermittivity^2; 
    reflectionCoefficient(:, reflectionOrderIndex) = [ ...
        (cosd(aor) - sqrt(B_v))/(cosd(aor) + sqrt(B_v)); ... % V
        (cosd(aor) - sqrt(B_h))/(cosd(aor) + sqrt(B_h))];    % H
end     
end

function incidentAngle = angleOfIncident(multipath)
diff25 = multipath(1,5:7)-multipath(1,2:4);
diff58 = multipath(1,5:7)-multipath(1,8:10);
diff25 = -diff25./norm(diff25);
diff58 = -diff58./norm(diff58);

dpAoI = dot(diff25,diff58);
incidentAngle(1) = 0.5*acosd(dpAoI);

if multipath(1,1) == 2
    diff85 = -diff58;
    diff811 = multipath(1,8:10)-multipath(1,11:13);
    diff811 = -diff811./norm(diff811);

    dpAoI = dot(diff85,diff811);
    incidentAngle(2) = 0.5*acosd(dpAoI);
end
end