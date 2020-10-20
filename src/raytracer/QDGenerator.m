function [output, mpcIdx, switchQd,paramsRotation] = QDGenerator(reflOrder,output,...
    arrayMaterials, number, materialLibrary, distance, freq, mpcIdx,...
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
    
    
    material = arrayMaterials(number,reflOrderId);
    if reflOrderId == reflOrder
        output(mpcIdx-1,21, 1:indStoc) = material;
    end
    
    pathLoss = PathlossQD(materialLibrary,...
        arrayMaterials(number,:),1);
    pathLossSum=pathLossSum+pathLoss;
end
output(mpcIdx-1,9, 1:indStoc) = output(mpcIdx-1,9,1:indStoc)-(pathLossSum);
pathGain = squeeze(output(mpcIdx-1,9,1:indStoc));

% i1 is for generating precursors and i2 is for generating postcursors.
for i1 = 1:2
    % In this step the material properties are extracted from material library
    % and different parameters are generated.
    
    if i1 == 1
        muk = materialLibrary.mu_k_Precursor(material);
        sigmak = materialLibrary.sigma_k_Precursor(material);
        muy = materialLibrary.mu_Y_Precursor(material);
        sigmay = materialLibrary.sigma_Y_Precursor(material);
        mul = materialLibrary.mu_lambda_Precursor(material);
        sigmal = materialLibrary.sigma_lambda_Precursor(material);
        muSigmas = materialLibrary.mu_SigmaS_Precursor(material);
        sigmaSigmas = materialLibrary.sigma_SigmaS_Precursor(material);
        n = Nprec;
        
    else
        muk = materialLibrary.mu_k_Postcursor(material);
        sigmak = materialLibrary.sigma_k_Postcursor(material);
        muy = materialLibrary.mu_Y_Postcursor(material);
        sigmay = materialLibrary.sigma_Y_Postcursor(material);
        mul = materialLibrary.mu_lambda_Postcursor(material);
        sigmal = materialLibrary.sigma_lambda_Postcursor(material);
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
        sigmaS = db2pow(normalRandomGenerator(muSigmas, sigmaSigmas,indStoc)); 
        kFactor = db2pow(normalRandomGenerator(muk, sigmak,indStoc));           
        gamma = normalRandomGenerator(muy, sigmay,indStoc);
        lambda = normalRandomGenerator(mul, sigmal,indStoc);
        sigmaAodEl = abs(normalRandomGenerator(muPhi, sigmaPhi,indStoc));
        sigmaAodAz = abs(normalRandomGenerator(muTheta, sigmaTheta,indStoc));
        sigmaAoaEl = abs(normalRandomGenerator(muPhi, sigmaPhi,indStoc));
        sigmaAoaAz = abs(normalRandomGenerator(muTheta, sigmaTheta,indStoc));
        
        lambda = abs(normalRandomGenerator(mul, sigmal,indStoc));
        tau_set = nan(n+1, indStoc);
        tau_set(1,:) = (distance/c) * 1e9; % convert to ns
        
        i = 1;
        
        % generates delays
        while i<=n
            
            diff = randomExponetialGenerator(lambda).';
            if i1 == 1
                tau_set(i+1,:) = tau_set(i,:)-diff;
            else
                tau_set(i+1,:) = tau_set(i,:)+diff;
            end
            output(mpcIdx+i-1,1,1:indStoc) = count-1;
            output(mpcIdx+i-1, 2:4,1:indStoc) = repmat(dod,1,1,indStoc);
            output(mpcIdx+i-1, 5:7,1:indStoc) = repmat(doa,1,1,indStoc);
            output(mpcIdx+i-1,8,1:indStoc) = tau_set(i+1,:)*1e-9;
            
            % output(count1,10)=180*atan(dod(2)/dod(1))/pi;
            % output(count1,11)=180*acos(dod(3)/norm(dod))/pi;
            % output(count1,12)=180*atan(doa(2)/doa(1))/pi;
            % output(count1,13)=180*acos(doa(3)/norm(doa))/pi;
            i=i+1;
        end
        
        pgCursor = db2pow(pathGain);
        % generates path loss
        
        for i=1:n
            
            if i~=0
                s = normalRandomGenerator(0,sigmaS);
                if i1==1
                       output(mpcIdx+i-1,9,1:indStoc) = pow2db((pgCursor./kFactor).*...
                           exp(((tau_set(i+1,:)-tau_set(1,:)).'./gamma)+s));
                else
                    output(mpcIdx+i-1,9,1:indStoc) =  pow2db((pgCursor./kFactor).*...
                           exp(-((tau_set(i+1,:)-tau_set(1,:)).'./gamma)+s));
                    if i == 1
                        output(indexReference,9,1:indStoc) =  20*log10(pgCursor-(pgCursor./kFactor));
                    end
                end
            end
        end
        
        mu = acosd(dod(3)/norm(dod));
        thetaCursor  = repmat(mu, 1, indStoc);
        
        % generates angular spread for Aod elevation
        Aod_el = zeros(n, indStoc);
        PG = zeros(n, indStoc);
        PG1 = zeros(n, indStoc);
        for i=1:n
            x = randomLaplaceGenerator(indStoc);
%             ran(i,:) = x;
            Aod_el(i,:) = mu+x;
            PG(i,:) = 10.^(output(mpcIdx+i-1,9,1:indStoc)/20);
            PG1(i,:) = output(mpcIdx+i-1,9,1:indStoc);
            %     output(count1+i-1,11)=randomLaplaceGenerator(mu,sigma_Aod_E);
        end
        PG_cursor = squeeze(10.^(output(count-1,9,1:indStoc)/20)).';
        PG_cursor1 = squeeze(output(count-1,9,1:indStoc)).';
        
        mu_1 = (sum(PG.*Aod_el)+(PG_cursor.*thetaCursor))...
            ./(sum(PG)+PG_cursor);
        sigma_1 = sqrt(((sum(PG.*((Aod_el-mu_1).*...
            (Aod_el-mu_1))))+(PG_cursor.*...
            ((thetaCursor-mu_1).^2)))./(sum(PG)+PG_cursor));
        s = sigmaAodEl.'./sigma_1;
        Aod_el1 = [thetaCursor;Aod_el];
        x = (s.^2)-((PG_cursor.*((thetaCursor-mu_1).^2))/...
            ((sum(PG.*((Aod_el-mu_1).*(Aod_el-mu_1))))+...
            (PG_cursor.*((thetaCursor-mu_1).^2))));
        a=((sum(PG.*((Aod_el-mu_1).*(Aod_el-mu_1))))+(PG_cursor.*((thetaCursor-mu_1).^2)));
        b=PG_cursor.*((thetaCursor-mu_1).^2);
        s1=abs(x.*(1+(b./a)));
        Aod_el=((Aod_el-mu_1).*sqrt(s1))+(mu_1);
        
        mu=squeeze(output(mpcIdx-1,10,1:indStoc));
        thetaCursor=mu.';
        % generates angular spread for Aod azimuth
        ran = [];
        for i=1:n
            x = randomLaplaceGenerator(indStoc);
            ran(i,:)=x;
            Aod_az(i,:)=mu+x;
            
        end
        
        mu_1=sum(PG.*Aod_az)./sum(PG);
        sigma_1=sqrt((sum(PG.*((Aod_az-mu_1).*(Aod_az-mu_1))))./sum(PG));
        s=sigmaAodAz.'./sigma_1;
        Aod_az1=[thetaCursor;Aod_az];
        x=(s.^2)-((PG_cursor.*((thetaCursor-mu_1).^2))./((sum(PG.*((Aod_az-mu_1).*(Aod_az-mu_1))))+(PG_cursor.*((thetaCursor-mu_1).^2))));
        a=((sum(PG.*((Aod_az-mu_1).*(Aod_az-mu_1))))+(PG_cursor.*((thetaCursor-mu_1).^2)));
        b=PG_cursor.*((thetaCursor-mu_1).^2);
        s1=abs(x.*(1+(b./a)));
        Aod_az=((Aod_az-mu_1).*sqrt(s1))+(mu_1);
        
        
        mu=squeeze(output(mpcIdx-1,12,1:indStoc));
        thetaCursor=mu.';
        % generates angular spread for Aoa azimuth
        
        ran = [];
        for i=1:n
            x = randomLaplaceGenerator(indStoc);
            ran(i,:)=x;
            Aoa_az(i,:)=mu+x;
%             theta_cursor=mu;
        end
        
        mu_1=sum(PG.*Aoa_az)./sum(PG);
        sigma_1=sqrt((sum(PG.*((Aoa_az-mu_1).*(Aoa_az-mu_1))))./sum(PG));
        s=sigmaAoaAz.'./sigma_1;
        Aoa_az1=[thetaCursor;Aoa_az];
        x=(s.^2)-((PG_cursor.*((thetaCursor-mu_1).^2))./((sum(PG.*((Aoa_az-mu_1).*(Aoa_az-mu_1))))+(PG_cursor.*((thetaCursor-mu_1).^2))));
        a=((sum(PG.*((Aoa_az-mu_1).*(Aoa_az-mu_1))))+(PG_cursor.*((thetaCursor-mu_1).^2)));
        b=PG_cursor.*((thetaCursor-mu_1).^2);
        s1=abs(x.*(1+(b./a)));
        Aoa_az=((Aoa_az-mu_1).*sqrt(s1))+(mu_1);
        
        mu=acosd(doa(3)/norm(doa));
        thetaCursor = repmat(mu,1,indStoc);

        ran = [];
        % generates angular spread for Aoa elevation
        for i=1:n
            x = randomLaplaceGenerator(indStoc);
            ran(i,:)=x;
            Aoa_el(i,:)=mu+x;
%             theta_cursor=mu;
        end
        
        mu_1=sum(PG.*Aoa_el)./sum(PG);
        sigma_1=sqrt((sum(PG.*((Aoa_el-mu_1).*(Aoa_el-mu_1))))./sum(PG));
        s=sigmaAoaEl.'./sigma_1;
        Aoa_el1=[thetaCursor;Aoa_el];
        x=(s.^2)-((PG_cursor.*((thetaCursor-mu_1).^2))/((sum(PG.*((Aoa_el-mu_1).*(Aoa_el-mu_1))))+(PG_cursor.*((thetaCursor-mu_1).^2))));
        a=((sum(PG.*((Aoa_el-mu_1).*(Aoa_el-mu_1))))+(PG_cursor.*((thetaCursor-mu_1).^2)));
        b=PG_cursor.*((thetaCursor-mu_1).^2);
        s1=abs(x.*(1+(b./a)));
        Aoa_el=((Aoa_el-mu_1).*sqrt(s1))+(mu_1);
        
        % sorting all QD parameters in 'QD' parametwr
        for i=1:n
            output(mpcIdx+i-1, 11, 1:indStoc)=Aod_el1(i+1,:);
            output(mpcIdx+i-1, 10, 1:indStoc)=Aod_az1(i+1,:);
            output(mpcIdx+i-1, 12, 1:indStoc)=Aoa_az1(i+1,:);
            output(mpcIdx+i-1, 13, 1:indStoc)=Aoa_el1(i+1,:);
            output(mpcIdx+i-1, 14:18, 1:indStoc)=output(mpcIdx-1, 14:18, 1:indStoc);
            
            vAngle_DoD = Aod_el1(i+1,:);
            hAngle_DoD = Aod_az1(i+1,:);
            dod_temp=reshape([sind(vAngle_DoD).*cosd(hAngle_DoD),sind(vAngle_DoD).*sind(hAngle_DoD),cosd(vAngle_DoD)], [], 3).';
                    
            dod_mpc_norm =  tau_set(i+1,:)/(tau_set(1,:))* norm(dod);
            dod_z =dod_mpc_norm.*cosd(vAngle_DoD);
            dod_y = sqrt(((dod_mpc_norm)^2-dod_z.^2)./(1./tand(hAngle_DoD).^2 +1));
            if hAngle_DoD>180
            dod_y = -dod_y;
            end
            dod_x = dod_y./tand(hAngle_DoD);
            dod_stochMPC = [dod_x, dod_y, dod_z];
            
            doa_mpc_norm =  tau_set(i+1,:)/(tau_set(1,:))* norm(doa);
            doa_z = doa_mpc_norm.*cosd(Aoa_el1(i+1,:));
            doa_y = sqrt(((doa_mpc_norm)^2-doa_z.^2)./(1./tand(Aoa_az1(i+1,:)).^2 +1));
            if Aoa_az1(i+1,:)>180
            doa_y = -doa_y;
            end
            doa_x = doa_y./tand(Aoa_az1(i+1,:));
            doa_stochMPC = [doa_x, doa_y, doa_z];
            
            
            vtx_along_dod = dot(repmat(vtx, indStoc,1).', -dod_temp);
            vrx_along_dod = dot(repmat(vTemp, indStoc,1).', -dod_temp);
            doppler_factor = freq * (vrx_along_dod - vtx_along_dod) / c;
            output(mpcIdx+i-1,20,1:indStoc) = doppler_factor;
            output(mpcIdx+i-1,18,1:indStoc) = rand*2*pi;
            output(mpcIdx+i-1,19,1:indStoc) = output(mpcIdx+i-1,9,1:indStoc).*(output(mpcIdx-1,19,1:indStoc)./output(mpcIdx-1,9,1:indStoc));
            paramsRotation(mpcIdx+i-1).dod = reshape(dod_stochMPC, [],3);
            paramsRotation(mpcIdx+i-1).doa = reshape(doa_stochMPC, [],3);
            paramsRotation(mpcIdx+i-1).TxVel = vtx_along_dod;
            paramsRotation(mpcIdx+i-1).RxVel = vrx_along_dod;
            switchQd = 1;
        end
        
        mpcIdx=mpcIdx+i;
        
    end
end

end
