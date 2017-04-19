clear;
clc;
close all;

L = 1; % how many unit cells in the horizontal direction
H = 5; % how many hexagons in the unit cell for a-gnr

lineatoms = 2 * H + 1; % half of the number of atoms in the unit cell

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
