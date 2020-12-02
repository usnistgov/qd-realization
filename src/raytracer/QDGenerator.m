function [output, mpcIdx, switchQd,paramsRotation] = QDGenerator(reflOrder,output,...
    arrayMaterials, indexMultipath, materialLibrary, distance, freq, mpcIdx,...
    dod, doa, vtx, vTemp, count, indexReference, paramsRotation, varargin)
%INPUT -
%order_of_R -  Order of reflection
%output - multipath parameters
%array_of_materials - material properties of every triangle as described in
%array of plane
%number - row number of array_of_materials that is being used
%Array_of_planes - Similar to Array of points. Each triangle occupies 4
%columns (plane equation). The first column has the order of reflection
%(o/p of treetraversal)
%Material_library - Similar to Array of points. Each triangle occupies 1
%triangle. The data is the row number of material from Material library
%array_of_materials - Similar to Array of points. Each triangle occupies 1
%triangle. The data is the row number of material from Material library
%distance is the total length of multipath
%freq - frequency of operation (center frequency of carrier wave)
%count1 - number of rows in output/ multipath
%dod - direction of departure vector
%doa - direction of arrival vector
%vtx is velocity of tx
%v_temp is relative velocity of rx wrt tx
%count - number of rows in output/ multipath
%
%OUTPUT -
%output - multipath parameters
%count - number of rows in output/ multipath
%switch_QD - whether QD output exists
%
% This part of code generates QD parameters of a select multipath


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
% Modified by: Mattia Lecci <leccimat@dei.unipd.it>, Improved access to MaterialLibrary
% Modified by: Neeraj Varshney <neeraj.varshney@nist.gov>, Included residual
% error in QD model

%% Varargin processing
p = inputParser;
addParameter(p,'indStoc',1)
addParameter(p,'Nprec', 3)
addParameter(p,'Npost', 16)
parse(p, varargin{:});
indStoc = p.Results.indStoc;
Nprec = p.Results.Nprec;
Npost = p.Results.Npost;

%% Init
c = getLightSpeed;
pathLossSum = 0;
switchQd = 0;

% if  switch_material==1 && QD_gen==1
% To demonstrate that even higher order reflections can be taken care of
% but for accurate results higer order reflections are not valid i.e.,
% above first order. For higher order reflections the physical measurements
% are in progress.
%  if  switch_material==1 && order_of_R==1
for reflOrderId = 1:reflOrder 
    material = arrayMaterials(indexMultipath,reflOrderId);

    if reflOrderId == reflOrder
        output(mpcIdx-1,21, 1:indStoc) = material;
    end
    
    pathLoss = PathlossQD(materialLibrary,...
        arrayMaterials(indexMultipath,:),1);
    pathLossSum=pathLossSum+pathLoss;
end
output(mpcIdx-1,9, 1:indStoc) = output(mpcIdx-1,9,1:indStoc)-(pathLossSum);
pathGain = squeeze(output(mpcIdx-1,9,1:indStoc));

% i1 is for generating precursors and i2 is for generating postcursors.
for cursorType = 1:2
    % In this step the material properties are extracted from material library
    % and different parameters are generated.
    
    if cursorType == 1
        muk = materialLibrary.mu_k_Precursor(material);
        sigmak = materialLibrary.sigma_k_Precursor(material);
        muy = materialLibrary.mu_Y_Precursor(material);
        sigmay = materialLibrary.sigma_Y_Precursor(material);
        muDelay = materialLibrary.mu_lambda_Precursor(material);
        sigmaDelay = materialLibrary.sigma_lambda_Precursor(material);
        muSigmas = materialLibrary.mu_SigmaS_Precursor(material);
        sigmaSigmas = materialLibrary.sigma_SigmaS_Precursor(material);
        n = Nprec;
        
    else
        muk = materialLibrary.mu_k_Postcursor(material);
        sigmak = materialLibrary.sigma_k_Postcursor(material);
        muy = materialLibrary.mu_Y_Postcursor(material);
        sigmay = materialLibrary.sigma_Y_Postcursor(material);
        muDelay = materialLibrary.mu_lambda_Postcursor(material);
        sigmaDelay = materialLibrary.sigma_lambda_Postcursor(material);
        muSigmas = materialLibrary.mu_SigmaS_Postcursor(material);
        sigmaSigmas = materialLibrary.sigma_SigmaS_Postcursor(material);
        n = Npost;
        
    end
    
    muTheta = materialLibrary.mu_sigmaThetaAZ(material);
    sigmaTheta = materialLibrary.sigma_sigmaThetaAZ(material);
    muPhi = materialLibrary.mu_sigmaThetaEL(material);
    sigmaPhi = materialLibrary.sigma_sigmaThetaEL(material);
    
    %%
    if  muk~=0
        sigmaS = db2pow(normalRandomGenerator(muSigmas, sigmaSigmas, indStoc)); 
        kFactor = db2pow(normalRandomGenerator(muk, sigmak,indStoc));           
        gamma = abs(normalRandomGenerator(muy, sigmay,indStoc));
        lambdaDelay = abs(normalRandomGenerator(muDelay, sigmaDelay,indStoc));
        sigmaAodEl = abs(normalRandomGenerator(muPhi, sigmaPhi,indStoc));
        sigmaAodAz = abs(normalRandomGenerator(muTheta, sigmaTheta,indStoc));
        sigmaAoaEl = abs(normalRandomGenerator(muPhi, sigmaPhi,indStoc));
        sigmaAoaAz = abs(normalRandomGenerator(muTheta, sigmaTheta,indStoc));
        tauSet = nan(n+1, indStoc);
        tauSet(1,:) = distance/c * 1e9; % convert to ns
        idxDiffuse = 1;
        
        %% Generates delays
        while idxDiffuse<=n            
            diff = randomExponetialGenerator(lambdaDelay).';

            if cursorType == 1
                tauSet(idxDiffuse+1,:) = tauSet(idxDiffuse,:)-diff;

            else
                tauSet(idxDiffuse+1,:) = tauSet(idxDiffuse,:)+diff;

            end

            output(mpcIdx+idxDiffuse-1,1,1:indStoc) = count-1;
            output(mpcIdx+idxDiffuse-1, 2:4,1:indStoc) = repmat(dod,1,1,indStoc);
            output(mpcIdx+idxDiffuse-1, 5:7,1:indStoc) = repmat(doa,1,1,indStoc);
            output(mpcIdx+idxDiffuse-1,8,1:indStoc) = tauSet(idxDiffuse+1,:)*1e-9;
            
            % output(count1,10)=180*atan(dod(2)/dod(1))/pi;
            % output(count1,11)=180*acos(dod(3)/norm(dod))/pi;
            % output(count1,12)=180*atan(doa(2)/doa(1))/pi;
            % output(count1,13)=180*acos(doa(3)/norm(doa))/pi;
            idxDiffuse=idxDiffuse+1;
        end
        
        %% generates path loss
        pgCursor = db2pow(pathGain);

        for idxDiffuse=1:n            
                s = normalRandomGenerator(0,sigmaS);

                if cursorType==1
                       output(mpcIdx+idxDiffuse-1,9,1:indStoc) = pow2db((pgCursor./kFactor).*...
                           exp(((tauSet(idxDiffuse+1,:)-tauSet(1,:)).'./gamma)+s));

                else
                    output(mpcIdx+idxDiffuse-1,9,1:indStoc) =  pow2db((pgCursor./kFactor).*...
                           exp(-((tauSet(idxDiffuse+1,:)-tauSet(1,:)).'./gamma)+s));

                    if idxDiffuse == 1 % probably to delete
                        output(indexReference,9,1:indStoc) =  20*log10(pgCursor-(pgCursor./kFactor));
                    end
                end
        end
        
        %% generates angular spread for Aod elevation
        
        mu = acosd(dod(3)/norm(dod));
        thetaCursor  = repmat(mu, 1, indStoc);
        
        aodEl = zeros(n, indStoc);
%         PG = zeros(n, indStoc);
%         PG1 = zeros(n, indStoc);
        for idxDiffuse=1:n
            x = randomLaplaceGenerator(indStoc);
%             ran(i,:) = x;
            aodEl(idxDiffuse,:) = mu+x;
%             PG(idxDiffuse,:) = 10.^(output(mpcIdx+idxDiffuse-1,9,1:indStoc)/20);
%             PG1(idxDiffuse,:) = output(mpcIdx+idxDiffuse-1,9,1:indStoc);
            %     output(count1+i-1,11)=randomLaplaceGenerator(mu,sigma_Aod_E);
        end
%         PG_cursor = squeeze(10.^(output(count-1,9,1:indStoc)/20)).';
%         PG_cursor1 = squeeze(output(count-1,9,1:indStoc)).';
        
%         mu_1 = (sum(PG.*aodEl)+(PG_cursor.*thetaCursor))...
%             ./(sum(PG)+PG_cursor);
%         sigma_1 = sqrt(((sum(PG.*((aodEl-mu_1).*...
%             (aodEl-mu_1))))+(PG_cursor.*...
%             ((thetaCursor-mu_1).^2)))./(sum(PG)+PG_cursor));
%         s = sigmaAodEl.'./sigma_1;
        aodEl = [thetaCursor;aodEl];
%         x = (s.^2)-((PG_cursor.*((thetaCursor-mu_1).^2))/...
%             ((sum(PG.*((aodEl-mu_1).*(aodEl-mu_1))))+...
%             (PG_cursor.*((thetaCursor-mu_1).^2))));
%         a=((sum(PG.*((aodEl-mu_1).*(aodEl-mu_1))))+(PG_cursor.*((thetaCursor-mu_1).^2)));
%         b=PG_cursor.*((thetaCursor-mu_1).^2);
%         s1=abs(x.*(1+(b./a)));
%         aodEl=((aodEl-mu_1).*sqrt(s1))+(mu_1);
        
        mu=squeeze(output(mpcIdx-1,10,1:indStoc));
        thetaCursor=mu.';
        aodAz = zeros(n, indStoc);

        % generates angular spread for Aod azimuth
%         ran = [];
        for idxDiffuse=1:n
            x = randomLaplaceGenerator(indStoc);
%             ran(idxDiffuse,:)=x;
            aodAz(idxDiffuse,:)=mu+x;             
        end
        
%         mu_1=sum(PG.*Aod_az)./sum(PG);
%         sigma_1=sqrt((sum(PG.*((Aod_az-mu_1).*(Aod_az-mu_1))))./sum(PG));
%         s=sigmaAodAz.'./sigma_1;
        aodAz=[thetaCursor;aodAz];
%         x=(s.^2)-((PG_cursor.*((thetaCursor-mu_1).^2))./((sum(PG.*((Aod_az-mu_1).*(Aod_az-mu_1))))+(PG_cursor.*((thetaCursor-mu_1).^2))));
%         a=((sum(PG.*((Aod_az-mu_1).*(Aod_az-mu_1))))+(PG_cursor.*((thetaCursor-mu_1).^2)));
%         b=PG_cursor.*((thetaCursor-mu_1).^2);
%         s1=abs(x.*(1+(b./a)));
%         Aod_az=((Aod_az-mu_1).*sqrt(s1))+(mu_1);
        
        
        mu=squeeze(output(mpcIdx-1,12,1:indStoc));
        thetaCursor=mu.';
        % generates angular spread for Aoa azimuth
        
%         ran = [];
        aoaAz = zeros(n, indStoc);

        for idxDiffuse=1:n
            x = randomLaplaceGenerator(indStoc);
%             ran(idxDiffuse,:)=x;
            aoaAz(idxDiffuse,:)=mu+x;
%             theta_cursor=mu;
        end
        
%         mu_1=sum(PG.*Aoa_az)./sum(PG);
%         sigma_1=sqrt((sum(PG.*((Aoa_az-mu_1).*(Aoa_az-mu_1))))./sum(PG));
%         s=sigmaAoaAz.'./sigma_1;
        aoaAz=[thetaCursor;aoaAz];
%         x=(s.^2)-((PG_cursor.*((thetaCursor-mu_1).^2))./((sum(PG.*((Aoa_az-mu_1).*(Aoa_az-mu_1))))+(PG_cursor.*((thetaCursor-mu_1).^2))));
%         a=((sum(PG.*((Aoa_az-mu_1).*(Aoa_az-mu_1))))+(PG_cursor.*((thetaCursor-mu_1).^2)));
%         b=PG_cursor.*((thetaCursor-mu_1).^2);
%         s1=abs(x.*(1+(b./a)));
%         Aoa_az=((Aoa_az-mu_1).*sqrt(s1))+(mu_1);
        
        mu=acosd(doa(3)/norm(doa));
        thetaCursor = repmat(mu,1,indStoc);
        aoaEl = zeros(n, indStoc);

%         ran = [];
        % generates angular spread for Aoa elevation
        for idxDiffuse=1:n
            x = randomLaplaceGenerator(indStoc);
%             ran(idxDiffuse,:)=x;
            aoaEl(idxDiffuse,:)=mu+x;
%             theta_cursor=mu;
        end
        
%         mu_1=sum(PG.*Aoa_el)./sum(PG);
%         sigma_1=sqrt((sum(PG.*((Aoa_el-mu_1).*(Aoa_el-mu_1))))./sum(PG));
%         s=sigmaAoaEl.'./sigma_1;
        aoaEl=[thetaCursor;aoaEl];
%         x=(s.^2)-((PG_cursor.*((thetaCursor-mu_1).^2))/((sum(PG.*((Aoa_el-mu_1).*(Aoa_el-mu_1))))+(PG_cursor.*((thetaCursor-mu_1).^2))));
%         a=((sum(PG.*((Aoa_el-mu_1).*(Aoa_el-mu_1))))+(PG_cursor.*((thetaCursor-mu_1).^2)));
%         b=PG_cursor.*((thetaCursor-mu_1).^2);
%         s1=abs(x.*(1+(b./a)));
%         Aoa_el=((Aoa_el-mu_1).*sqrt(s1))+(mu_1);
        
        % sorting all QD parameters in 'QD' parametwr
        for idxDiffuse=1:n
            output(mpcIdx+idxDiffuse-1, 11, 1:indStoc)=aodEl(idxDiffuse+1,:);
            output(mpcIdx+idxDiffuse-1, 10, 1:indStoc)=aodAz(idxDiffuse+1,:);
            output(mpcIdx+idxDiffuse-1, 12, 1:indStoc)=aoaAz(idxDiffuse+1,:);
            output(mpcIdx+idxDiffuse-1, 13, 1:indStoc)=aoaEl(idxDiffuse+1,:);
            output(mpcIdx+idxDiffuse-1, 14:18, 1:indStoc)=output(mpcIdx-1, 14:18, 1:indStoc);
            
            vAngleDoD = aodEl(idxDiffuse+1,:);
            hAngleDoD = aodAz(idxDiffuse+1,:);
            dod_temp=reshape([sind(vAngleDoD).*cosd(hAngleDoD),sind(vAngleDoD).*sind(hAngleDoD),cosd(vAngleDoD)], [], 3).';
                    
            dod_mpc_norm =  tauSet(idxDiffuse+1,:)/(tauSet(1,:))* norm(dod);
            dod_z =dod_mpc_norm.*cosd(vAngleDoD);
            dod_y = sqrt(((dod_mpc_norm)^2-dod_z.^2)./(1./tand(hAngleDoD).^2 +1));
            if hAngleDoD>180
            dod_y = -dod_y;
            end
            dod_x = dod_y./tand(hAngleDoD);
            dod_stochMPC = [dod_x, dod_y, dod_z];
            
            doa_mpc_norm =  tauSet(idxDiffuse+1,:)/(tauSet(1,:))* norm(doa);
            doa_z = doa_mpc_norm.*cosd(aoaEl(idxDiffuse+1,:));
            doa_y = sqrt(((doa_mpc_norm)^2-doa_z.^2)./(1./tand(aoaAz(idxDiffuse+1,:)).^2 +1));
            if aoaAz(idxDiffuse+1,:)>180
            doa_y = -doa_y;
            end
            doa_x = doa_y./tand(aoaAz(idxDiffuse+1,:));
            doa_stochMPC = [doa_x, doa_y, doa_z];
            
            
            vtx_along_dod = dot(repmat(vtx, indStoc,1).', -dod_temp);
            vrx_along_dod = dot(repmat(vTemp, indStoc,1).', -dod_temp);
            doppler_factor = freq * (vrx_along_dod - vtx_along_dod) / c;
            output(mpcIdx+idxDiffuse-1,20,1:indStoc) = doppler_factor;
            output(mpcIdx+idxDiffuse-1,18,1:indStoc) = rand*2*pi;
            output(mpcIdx+idxDiffuse-1,19,1:indStoc) = output(mpcIdx+idxDiffuse-1,9,1:indStoc).*(output(mpcIdx-1,19,1:indStoc)./output(mpcIdx-1,9,1:indStoc));
            paramsRotation(mpcIdx+idxDiffuse-1).dod = reshape(dod_stochMPC, [],3);
            paramsRotation(mpcIdx+idxDiffuse-1).doa = reshape(doa_stochMPC, [],3);
            paramsRotation(mpcIdx+idxDiffuse-1).TxVel = vtx_along_dod;
            paramsRotation(mpcIdx+idxDiffuse-1).RxVel = vrx_along_dod;
            switchQd = 1;
        end
        
        mpcIdx=mpcIdx+idxDiffuse;
        
    end
end

end
