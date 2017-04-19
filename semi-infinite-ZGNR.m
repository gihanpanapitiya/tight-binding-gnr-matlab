
% tight binding calculations for zgnr

clear;
clc;
close all;

unitcellatms = 10;

L = 8;

y(1,1) = 0;
x(1,1) = 1.23;
xswitch = 0;

for i = 2:unitcellatms
  
    if mod(i,2)==1
        y(1,i) = y(1,i-1) + 1.42;
    else
        y(1,i) = y(1,i-1) + 0.71;
      
        x(1,i) = xswitch;
      
        if (i+1)<= unitcellatms
            x(1,i+1) = xswitch;
        end
      
        if xswitch == 0
            xswitch =1.23;
        else
            xswitch = 0;
        end
    end
  
end

for i = 2: L
   for j = 1: unitcellatms
       x(i,j) = x(i-1,j) + 2.46;
       y(i,j) = y(1,j);
      
   end
end


figure;
for i = 1:L
    plot (x(i,:),y(i,:), 'bo','MarkerFaceColor','b');
    hold on;
end
axis equal;


% setting up the nearest neighbours

nn(1,1) = 2;
nn(1,2) = 2;
nn(1,3) = 0;

nn(unitcellatms,1) = unitcellatms-1;
nn(unitcellatms,2) = unitcellatms-1;
nn(unitcellatms,3) = 0;

for i = 2:unitcellatms-1
  
    if mod(i,2) ==0
        nn(i,1) = i-1;
        nn(i,2) = i-1;
        nn(i,3) = i+1;
      
    else
        nn(i,1) = i+1;
        nn(i,2) = i+1;
        nn(i,3) = i-1;
      
    end
end

X = x(1,:); Y = y(1,:);

a = 2.46; 
e2p =0; t=-1; s=.2;
kx = -1*pi/a:.01:1*pi/a;
ky = -.05*pi/a:.01:.05*pi/a;

a1 = 1.23; a2 = 0.71; a3 = 1.42;

for m = 1: length(kx)
    for n = 1:length(ky)
      
        KX = kx(m);
        KY = ky(n);
      
        for i = 1:unitcellatms
          
            for j = 1:2:3
              
                nene = nn(i,j);

                if nene ~= 0
                    num = (Y(nene)-Y(i));
                    denom = (X(nene)-X(i));
          
                if (denom > 0 && num > 0)
                    H(i,nene) = exp(1i*KX*a1 + 1i*KY*a2)+exp(-1i*KX*a1 + 1i*KY*a2);
                    S(i,nene) = s * H(i,nene);

                elseif (denom < 0 && num > 0)
                    H(i,nene) = exp(-1i*KX*a1 + 1i*KY*a2)+exp(1i*KX*a1 + 1i*KY*a2);
                    S(i,nene) = s * H(i,nene);
                  
                elseif (denom < 0 && num < 0)
                    H(i,nene) = exp(-1i*KX*a1 - 1i*KY*a2)+exp(1i*KX*a1 - 1i*KY*a2);
                    S(i,nene) = s * H(i,nene);
                  
                elseif (num <0 && denom ==0)
                    H(i,nene) = exp(-1i*KY*a3);
                    S(i,nene) = s * H(i,nene);
                  
                elseif (denom > 0 && num < 0)
                    H(i,nene) = exp(1i*KX*a1 - 1i*KY*a2) + exp(-1i*KX*a1 - 1i*KY*a2) ;
                    S(i,nene) = s * H(i,nene);
                  
                elseif (num>0 && denom ==0)
                    H(i,nene) = exp(1i*KY*a3);
                    S(i,nene) = s * H(i,nene);
                  
                end
              
                end
              
            end
      
            H(i,i) = 0;
            S(i,i) = 1;
          
        end
      
        ee = eig(t*H,S);


        for k = 1:unitcellatms
            E(m,n,k) = ee(k);
        end
      
    end
end
%
figure;
for k = 1 : unitcellatms
  
    plot(kx, E(:,:,k));
    hold on;
  
end
