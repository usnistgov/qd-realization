function angleVector = angleToVector(az, el, delay)

angleVectorNorm = getLightSpeed*delay;


angleVectorZ =angleVectorNorm.*cosd(el);
signY = double(az<180);
signY(signY == 0 ) = -1;
angleVectorY =signY.*sqrt((angleVectorNorm.^2-angleVectorZ.^2)./(1./tand(az).^2 +1));
angleVectorX = angleVectorY./tand(az);
angleVector = [angleVectorX, angleVectorY, angleVectorZ];


