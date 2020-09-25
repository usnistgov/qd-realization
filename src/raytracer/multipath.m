function [QD, switchQD, output, multipath, indexMultipath, indexQD,varargout] =...
    multipath(ArrayOfPlanes, ArrayOfPoints, Rx, Tx, CADOutput,...
    numberOfRowsArraysOfPlanes, MaterialLibrary, arrayOfMaterials,...
    switchMaterial, velocityTx, velocityRx, PolarizationSwitch,...
    PolarizationTx, AntennaOrientationTx, PolarizationRx,...
    AntennaOrientationRx, switchCrossPolarization, QDGeneratorSwitch, frequency, varargin)
%INPUT -
%ArrayOfPoints - combinations of multiple triangles, every row is a unique
%combination. every triangle occupies 9 columns (3 vertices). (o/p of
%treetraversal)
%ArrayOfPlanes - Similar to Array of points. Each triangle occupies 4
%columns (plane equation). The first column has the order of reflection
%(o/p of treetraversal)
%Rx - Rx position
%Tx - Tx position
%CADop - CAD output
%number1 -
%MaterialLibrary - Similar to Array of points. Each triangle occupies 1
%triangle. The data is the row number of material from Material library
%arrayOfMaterials - Similar to Array of points. Each triangle occupies 1
%triangle. The data is the row number of material from Material library
%switchMaterial - whether triangle materials properties are present
% vtx, vrx are velocities of tx and rx respectively
% PolarizationSwitch - switch to enable or disable polarization module
% PolarizationTx/ PolarizationRx - Tx/Rx Polarization
% AntennaOrientationTx/ AntennaOrientationRx - Tx/Rx antenna
% oreientation
%switchCrossPolarization - a boolean to describe whether cross polarization is selected
%or not. 1 means there is cross polarization and 0 means there is no cross
%polarization
% QDGeneratorSwitch - Switch to turn ON or OFF the Qausi dterministic module
% 1 = ON, 0 = OFF
% frequency: the carrier frequency at which the system operates
%
%OUTPUT -
%QD - output to be plotted on f2 channel plot
%switchQD - whether QD output exists
%output - multipath parameters
%multipath - output to be plottd on f1 multipath plot
%indexMultipath - number of rows in multipath
%countQD - number of QD component
%
% The phase information in case of presence of polarization information and is
% encoded in the Jones vector. In case of absence of polarization, order of
% reflection is multiplied by pi to give phase shift


% -------------Software Disclaimer---------------
%
% NIST-developed software is provided by NIST as a public service. You may use, copy
% and distribute copies of the software in any medium, provided that you keep intact this
% entire notice. You may improve, modify and create derivative works of the software or
% any portion of the software, and you may copy and distribute such modifications or
% works. Modified works should carry a notice stating that you changed the software
% and should note the date and nature of any such change. Please explicitly
% acknowledge the National Institute of Standards and Technology as the source of the
% software.
%
% NIST-developed software is expressly provided "AS IS." NIST MAKES NO
% WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT OR ARISING BY
% OPERATION OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
% WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
% NON-INFRINGEMENT AND DATA ACCURACY. NIST NEITHER REPRESENTS
% NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE
% UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE
% CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS
% REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF,
% INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY,
% RELIABILITY, OR USEFULNESS OF THE SOFTWARE.
%
% You are solely responsible for determining the appropriateness of using
% and distributing the software and you assume all risks associated with its use, including
% but not limited to the risks and costs of program errors, compliance with applicable
% laws, damage to or loss of data, programs or equipment, and the unavailability or
% interruption of operation. This software is not intended to be used in any situation
% where a failure could cause risk of injury or damage to property. The software
% developed by NIST employees is not subject to copyright protection within the United
% States.
%
% Modified by: Mattia Lecci <leccimat@dei.unipd.it>, Used MATLAB functions instead of custom ones,
%    vectorized code, improved access to MaterialLibrary

var_struct = {'indStoc', 'qTx', 'qRx'}; 
for k = 1:2:length(varargin)
    if (~any(strcmp(varargin{k}, var_struct)))
        warning(['Cannot specify "', varargin{k}, '" as input value - it will be discarted']);
    end
    eval([varargin{k},' = varargin{k+1};'])
end
if ~exist('indStoc','var'),     indStoc=1;        end
if ~exist('qTx','var'),     qTx.cTx = Tx;  qTx.euc = zeros(1,3);   end
if ~exist('qRx','var'),     qRx.cRx = Rx;  qRx.euc = zeros(1,3);   end

switchQD = 0;
QD = [];
indexMultipath = 1;
indexOutput = 1;
indexQD = 1;
sizeArrayOfPlanes = size(ArrayOfPlanes);
paramsRotation = [];

output = zeros(sizeArrayOfPlanes(1),21,indStoc);
multipath = [];
LIGHTVELOCITY = 3e8;
wavelength = LIGHTVELOCITY / frequency;
if numberOfRowsArraysOfPlanes>0
    orderOfReflection = ArrayOfPlanes(1,1);
    
    % the for loop iterates through all the rows of ArrayOfPlanes,
    % ArrayOfPoints and provides a single row as input to
    % singleMultipathGenerator function
    % QD model is present in  this loop
    multipath = zeros(numberOfRowsArraysOfPlanes,orderOfReflection * 3 + 1);
    for iterateNumberOfRowsArraysOfPlanes = 1:numberOfRowsArraysOfPlanes
        
        
        indexOrderOfReflection = 1;
        multipath(indexMultipath, (indexOrderOfReflection-1)*3 + 1) = orderOfReflection;
        multipath(indexMultipath, (indexOrderOfReflection-1)*3 + 1 + (1:3)) = Rx;
        Reflected = Rx;
        
        % a single row of ArrayOfPlanes,ArrayOfPoints is fed to
        % singleMultipathGenerator function to know whether a path exists. If a
        % path exists then what are vectors that form the path (stored in
        % multipath parameter)
        
        PolarizationSwitchTemporary = 1;
        if PolarizationSwitch == 1 && switchMaterial == 1
            for orderOfReflectionTemporary = 1:orderOfReflection
                dielectricConstant = MaterialLibrary.DielectricConstant(...
                    arrayOfMaterials(indexMultipath,orderOfReflectionTemporary));
                if dielectricConstant ~= 0
                    nt_array(orderOfReflectionTemporary) = dielectricConstant;
                else
                    nt_array = [];
                    PolarizationSwitchTemporary = 0;
                    break
                end
            end
        else
            nt_array = [];
            PolarizationSwitchTemporary = 0;
        end
        
        [isMPC,~,dodNoRot,doaNoRot,multipath,distance,dopplerFactor,PathLoss,...
            ~,~,~,~,velocityTemp] = singleMultipathGenerator...
            (iterateNumberOfRowsArraysOfPlanes,orderOfReflection,indexOrderOfReflection,ArrayOfPlanes,...
            ArrayOfPoints,Reflected,Rx,Tx,CADOutput,...
            multipath,indexMultipath,velocityTx,velocityRx,PolarizationSwitchTemporary,...
            PolarizationTx,AntennaOrientationTx,PolarizationRx,...
            AntennaOrientationRx,nt_array,switchCrossPolarization);
            dod=  coordinateRotation(dodNoRot,[0 0 0], qTx.euc, 'frame');
            doa=  coordinateRotation(doaNoRot,[0 0 0], qRx.euc, 'frame'); 
        if isMPC == 1
            for i = 1:indexMultipath - 1
                switch3 = 1;
                for j = 1:(orderOfReflection * 3) + 6
                    switch3 = switch3 && (multipath(i,j) == multipath(indexMultipath,j));
                end
                isMPC = isMPC && (~switch3);
            end
        end
        
        % the delay, AoA, AoD, path loss of the path are stored in output parameter
        if  isMPC == 1
            
            output(indexOutput,1,1:indStoc) = indexMultipath;
            % dod - direction of departure
            output(indexOutput,2:4,1:indStoc) = repmat(dod,1,1,indStoc);
            % doa - direction of arrival
            output(indexOutput,5:7,1:indStoc) = repmat(doa,1,1,indStoc);
            % Time delay
            output(indexOutput,8,1:indStoc) = distance / LIGHTVELOCITY;
            % Friis transmission loss
            if PathLoss(1) < 0
                output(indexOutput,9,1:indStoc) = 20*log10(wavelength / (4*pi*distance)) + PathLoss(1);
            else
                output(indexOutput,9,1:indStoc) = 20*log10(wavelength / (4*pi*distance)) - 10;
            end
            % Aod azimuth
            output(indexOutput,10,1:indStoc) = mod(atan2d(dod(2),dod(1)), 360);
            % Aod elevation
            output(indexOutput,11,1:indStoc) = acosd(dod(3) / norm(dod));
            % Aoa azimuth
            output(indexOutput,12,1:indStoc) = mod(atan2d(doa(2),doa(1)), 360);
            % Aoa elevation
            output(indexOutput,13,1:indStoc) = acosd(doa(3) / norm(doa));
            output(indexOutput,18,1:indStoc) = orderOfReflection*pi;% + dopplerFactor*delay;
            output(indexOutput,20,1:indStoc) = dopplerFactor * frequency;
            paramsRotation(indexOutput).dod = dod;
            paramsRotation(indexOutput).doa = doa;
            paramsRotation(indexOutput).TxVel = velocityTx;
            paramsRotation(indexOutput).RxVel = velocityRx;
            indexReference = indexOutput;
            indexMultipath = indexMultipath + 1;
            indexOutput = indexOutput + 1;
            output(indexOutput - 1,21,1:indStoc) = 0;
            
            % refer to "multipath - WCL17_revised.pdf" in this folder for QD model
            if  switchMaterial == 1 && QDGeneratorSwitch == 1
                [output,indexOutput,switchQD,paramsRotation] = QDGenerator(orderOfReflection,...
                    output,arrayOfMaterials,iterateNumberOfRowsArraysOfPlanes,MaterialLibrary,distance,...
                    frequency,indexOutput,dod,doa,velocityTx,velocityTemp,indexMultipath,indexReference,paramsRotation, 'indStoc', indStoc);
            end
        end
        
    end
    
    switch2 = 1;
    
    try
        ioi = output(indexMultipath,1,1:indStoc)==0;
    catch
        switch2 = 0;
        indexMultipath = indexMultipath -1;
    end
    
    if switch2==1
        if output(indexMultipath,1,1:indStoc)==0
            indexMultipath = indexMultipath-1 ;
        end
    end
    
    indexQD = indexOutput - 1;
    output1 = output;
    output = output1(1:indexQD,1:21,1:indStoc);
    varargout{1} = paramsRotation;

    mp1 = multipath;
    multipath = [];
    
    if indexMultipath>=1
        multipath = mp1(1:indexMultipath,:);
    end
end

end