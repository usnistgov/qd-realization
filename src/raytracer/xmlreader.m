function [CADOutput, materialSwitch] = xmlreader(filename,...
    MaterialLibrary, referencePoint, r, IndoorSwitch)
%XML redaer extracts the information of CAD file (AMF). The input of the
%function is filename, the material database with all the material
%parameters (Material_library), reference point (referencePoint) and distance limitation(r)
%The output is extracted triangles (CADop),
%and a boolean to know whether the material information is
%present in the CAD file (switch1)


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
% anddistributing the software and you assume all risks associated with its use, including
% but not limited to the risks and costs of program errors, compliance with applicable
% laws, damage to or loss of data, programs or equipment, and the unavailability or
% interruption of operation. This software is not intended to be used in any situation
% where a failure could cause risk of injury or damage to property. The software
% developed by NIST employees is not subject to copyright protection within the United
% States.
%
% Modified by: Mattia Lecci <leccimat@dei.unipd.it>, Used MATLAB functions instead of custom ones,
%    improved MaterialLibrary access, readibility, performance in general
% Modified by: Neeraj Varshney <neeraj.varshney@nist.gov>, support multiple
% objects and different length units in amf file


s = xml2struct(filename);

% Probing whether material information is present or not
if isfield(s.amf,'material')
    materialSwitch = 1;
    
    sizeMaterials1 = size(s.amf.material');
    if sizeMaterials1(2)>1 &&  sizeMaterials1(1)==1
        sizeMaterials = sizeMaterials1;
    else
        sizeMaterials = sizeMaterials1(1);
    end
    
else
    materialSwitch = 0;
    
end
%% Iterating through all the subdivisions (volumes) and extracting the triangle information

% countRows=1; 
CADOutput = [];
objLength=length(s.amf.object);                                                            
for iobjLength = 1:objLength                            % For multiple objects          
    if objLength == 1
        volume=s.amf.object.mesh.volume';                                     
        sizeVolume=size(volume);
        stemp = s.amf.object;                                                 
    else
        volume=s.amf.object{1,iobjLength}.mesh.volume';                                     
        sizeVolume=size(volume);
        stemp = s.amf.object{1,iobjLength};                                                 
    end
    for iterateVolume=1:sizeVolume
        if sizeVolume(1)~=1
            triangles=stemp.mesh.volume{1,iterateVolume}.triangle';                     
        else
            triangles=stemp.mesh.volume.triangle';                                      
        end

        if materialSwitch==1
            if sizeVolume(1)~=1
                materialid=stemp.mesh.volume{1,iterateVolume}.Attributes.materialid;   
            else
                materialid=stemp.mesh.volume.Attributes.materialid;                     
            end

            for iterateMaterials=1:sizeMaterials
                if sizeMaterials~=1
                    if str2double(materialid)==str2double(s.amf.material{1,iterateMaterials}.Attributes.id)
                        Material=s.amf.material{1,iterateMaterials}.metadata.Text;
                    end

                elseif sizeVolume(1)==1 && sizeMaterials == 1
                    if str2double(materialid)==str2double(s.amf.material.Attributes.id)
                        Material=s.amf.material.metadata.Text;
                    end

                end
            end
        end
        %% Extracting the vertices information of the triangles

        sizeTriangle=size(triangles);
        CADOutputTemp = [];
        for iterateTriangles=1:sizeTriangle(1)
%             indexCADOutput=countRows;
            if isfield(s.amf,'Attributes')
                switch s.amf.Attributes.unit
                    case 'millimeter'
                        v1 = getTriangleVertex(stemp,iterateVolume,iterateTriangles,'v1',sizeVolume)*0.001;   
                        v2 = getTriangleVertex(stemp,iterateVolume,iterateTriangles,'v2',sizeVolume)*0.001;   
                        v3 = getTriangleVertex(stemp,iterateVolume,iterateTriangles,'v3',sizeVolume)*0.001;   
                    case 'inch'
                        v1 = getTriangleVertex(stemp,iterateVolume,iterateTriangles,'v1',sizeVolume)*0.0254;   
                        v2 = getTriangleVertex(stemp,iterateVolume,iterateTriangles,'v2',sizeVolume)*0.0254;   
                        v3 = getTriangleVertex(stemp,iterateVolume,iterateTriangles,'v3',sizeVolume)*0.0254;   
                    case 'foot'
                        v1 = getTriangleVertex(stemp,iterateVolume,iterateTriangles,'v1',sizeVolume)*0.3048;   
                        v2 = getTriangleVertex(stemp,iterateVolume,iterateTriangles,'v2',sizeVolume)*0.3048;   
                        v3 = getTriangleVertex(stemp,iterateVolume,iterateTriangles,'v3',sizeVolume)*0.3048;   
                    case 'micrometer'
                        v1 = getTriangleVertex(stemp,iterateVolume,iterateTriangles,'v1',sizeVolume)*1e-6;   
                        v2 = getTriangleVertex(stemp,iterateVolume,iterateTriangles,'v2',sizeVolume)*1e-6;   
                        v3 = getTriangleVertex(stemp,iterateVolume,iterateTriangles,'v3',sizeVolume)*1e-6;   
                    case 'meter'
                        v1 = getTriangleVertex(stemp,iterateVolume,iterateTriangles,'v1',sizeVolume)*1;   
                        v2 = getTriangleVertex(stemp,iterateVolume,iterateTriangles,'v2',sizeVolume)*1;   
                        v3 = getTriangleVertex(stemp,iterateVolume,iterateTriangles,'v3',sizeVolume)*1;                   
                    otherwise
                        error('xmlreader does not support this unit.');
                end
            else
               error('Length unit is missing in the xml/amf file. Add <amf  unit="?" in Line 1>');
            end
            % Calculating the plane equation of triangles

            vector1 = v2 - v3;
            vector2 = -(v2 - v1);
            % Multiply with -1 for 'example.xml', 'sphere.xml','material_prism.xml'
            normal = cross(vector2,vector1) * (1-(2*IndoorSwitch));
            normal = round(normal/norm(normal),4);
            vector3 = v2;
            % for box. remove for others
            D = -dot(normal,vector3);

            % Storing Material information in output if the material exists in the material database
            if materialSwitch==1
                switch2=0;
                for iterateMaterials=1:size(MaterialLibrary,1)
                    if strcmpi(MaterialLibrary.Reflector{iterateMaterials},Material)
                        CADOutputTemp(14) = MaterialLibrary.PrimaryKey(iterateMaterials);
                        switch2=1;
                    end
                end

                % Storing triangle vertices and plane equations in output
                % Part where output file is created. It contains the triangle vertices
                % in first nine columns, plane equations in the next four columns
                if switch2==0
                    materialSwitch=0;
                end

            end

            CADOutputTemp(1:3) = round(v1,6);
            CADOutputTemp(4:6) = round(v2,6);
            CADOutputTemp(7:9) = round(v3,6);
            CADOutputTemp(10:12) = round(normal,4);
            CADOutputTemp(13) = round(D,4);

            % We are using distance limitation at this step
            if r==0
                [switchDistance] = 1;
            else
                [switchDistance] = verifydistance(r,referencePoint,CADOutputTemp,1);
            end

            % If the triangles are within the given distance we increase the count,
            % else the next triangle will replace the present row (as count remains constant)
            if switchDistance==1
%                 countRows=countRows+1;
               CADOutput = [CADOutput;CADOutputTemp];
            end

        end
    end
end
end


%% Utils
function v = getTriangleVertex(stemp,volumeIdx,triangIdx,vertexIdx,sizeVolume)                      

if sizeVolume(1)~=1
    vertex = str2double(stemp.mesh.volume{1,volumeIdx}.triangle{1,triangIdx}.(vertexIdx).Text)+1;   
else
    vertex = str2double(stemp.mesh.volume.triangle{1,triangIdx}.(vertexIdx).Text)+1;                
end

x = str2double(stemp.mesh.vertices.vertex{1,vertex}.coordinates.x.Text);                            
y = str2double(stemp.mesh.vertices.vertex{1,vertex}.coordinates.y.Text);                     
z = str2double(stemp.mesh.vertices.vertex{1,vertex}.coordinates.z.Text);                     
v = [x,y,z];

end