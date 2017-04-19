 % tight binding calculations for agnr
 clear;
 clc;
 close all;

 L = 2;
 H = 4;

 lineatoms = 2 * H + 1;

 for i = 1:lineatoms

     if mod(i,2)==1
         x(1,i) = .71;
     else
         x(1,i) = 0;
     end

      y(1,i) = (i-1)*1.23;
 end

 index = 1;

 for i = lineatoms +1 : 2*lineatoms

     if mod(i,2)==0
         x(1,i) = 2.13;
     else
         x(1,i)= 2.84;
     end

     y(1,i) = y(1,i -index);

     index = index+2;


 end



 for i = 2:L
     for j = 1: 2*lineatoms
     x(i,j) = x(1,j) + 4.26*(i-1);
     y(i,j) = y(1,j);
     end
 end

 figure;
 for i = 1:L
     plot (x(i,:),y(i,:), 'bo','MarkerFaceColor','b');
     hold on;
 end
 axis equal;

 z = size(x);
 lastatm = z(2);
 atomnums = 1: lastatm;
 corneratms = [1 lastatm/2 lastatm/2+1 lastatm];

 nn(1,1) = 1+1;
 nn(1,2) = lastatm;
 nn(lastatm,1) = 1;
 nn(lastatm,2) = lastatm-1;


 nn(corneratms(2),1) = corneratms(2)-1;
 nn(corneratms(2),2) = corneratms(2)+1;
 nn(corneratms(3),1) = corneratms(3)-1;
 nn(corneratms(3),2) = corneratms(3)+1;

 X = x(1,:); Y = y(1,:);
  
 batm = lastatm-3;

 for i = 2:lastatm/2 -1
    
     nn(i,1) = i-1;
     nn(i,2) = i+1;
     nn(i,3) = i+batm;

     nnf(i,1) = i-1;
     nnf(i,2) = i+1;
      
    batm = batm -2;
end

batm = 3;
  
for i = lastatm/2+2:lastatm-1
      
    nn(i,1) = i-batm;
    nn(i,2) = i-1;
    nn(i,3) = i+1;
  
    batm = batm + 2;
end
  
a = 2.46; 
e2p =0; t=-1; s=.2;
kx = -1*pi/a:.01:1*pi/a;
ky = -.05*pi/a:.01:.05*pi/a;

a1 = 0.71; a2 = 1.23; a3 = 1.42;

for m = 1: length(kx)
    for n = 1:length(ky)
      
        KX = kx(m);
        KY = ky(n);
      
        for i = 1:2*lineatoms
  
            for j = 1:3
                              
                nene = nn(i,j); % three nearest neighbours
          
                if nene ~= 0
                    num = (Y(nene)-Y(i));
                    denom = (X(nene)-X(i));
          
                if (denom > 0 && num > 0)
                    H(i,nene) = exp(1i*KX*a1 + 1i*KY*a2);
                    S(i,nene) = s * H(i,nene);
                elseif (denom < 0 && num < 0)
                    H(i,nene) = exp(-1i*KX*a1 - 1i*KY*a2);
                    S(i,nene) = s * H(i,nene);
                elseif (denom < 0 && num > 0)
                    H(i,nene) = exp(-1i*KX*a1 + 1i*KY*a2);
                    S(i,nene) = s * H(i,nene);
                elseif (denom > 0 && num < 0)
                    H(i,nene) = exp(1i*KX*a1 - 1i*KY*a2);
                    S(i,nene) = s * H(i,nene);
                elseif (denom > 0 && num == 0 && denom ==1.42)
                    H(i,nene) = exp(1i*KX*a3);
                    S(i,nene) = s * H(i,nene);
                elseif (denom > 0 && num == 0 && denom == 2*1.42)
                    H(i,nene) = exp(-1i*KX*a3);
                    S(i,nene) = s * H(i,nene);
                elseif (denom < 0 && num == 0 && denom == -1.42)
                    H(i,nene) = exp(-1i*KX*a3);
                    S(i,nene) = s * H(i,nene);
                elseif (denom < 0 && num == 0 && denom == -2*1.42)
                    H(i,nene) = exp(1i*KX*a3);
                    S(i,nene) = s * H(i,nene);
                end
              
                end
            end
            H(i,i) = 0;
            S(i,i) = 1;
        end
      
        ee = eig(t*H,S); % Eigen value calculation
        % ee = eig( H );


        for k = 1:2*lineatoms
            E(m,n,k) = ee(k);
        end
      
    end
end

figure;
for k = 1 :2*lineatoms
  
    plot(kx, E(:,:,k));
    hold on;
  
end
