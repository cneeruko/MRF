close all
clear all

w1=2;
w2=4;
h1=3;
h2=2;
w3=2;
h3=2;

% position_A=[2*rand(1) 0];
% bx=10*rand(1);
% while ((bx<position_A(1)-w2) ||(bx>position_A(1)+w1))
%    bx=10*rand(1);
% end
% position_B=[bx position_A(2)+h1];
% 
% position_C=[3*rand(1)+ position_A(1)+w1 0];

position_A=[1,0];
position_B=[0,3];
position_C=[5,0];



figure(1)
A=rectangle('Position',[position_A(1) position_A(2) w1 h1],'EdgeColor','b','LineWidth',3);
hold on 
B=rectangle('Position',[position_B(1) position_B(2) w2 h2],'EdgeColor','r','LineWidth',3);
hold on 
C=rectangle('Position',[position_C(1) position_C(2) w3 h3],'EdgeColor','g','LineWidth',3);
axis equal
axis([-10 10 0 10])


%center_A=[2,1.5]
%center_B=[2,4]


positions_for_A(1,:)=[position_A(1) normrnd(position_A(1),0.5,1,4)];
positions_for_B(1,:)=[position_B(1) normrnd(position_B(1),0.5,1,4)];
positions_for_A(2,:)=[position_A(2) normrnd(position_A(2),0.5,1,4)];
positions_for_B(2,:)=[position_B(2) normrnd(position_B(2),0.5,1,4)];
positions_for_C(1,:)=[position_C(1) normrnd(position_C(1),0.5,1,4)];
positions_for_C(2,:)=[position_C(2) normrnd(position_C(2),0.5,1,4)];

% positions_for_B= [0,0.600142508808630,0.593930795648753,0.249517740562795,1.80886262051925,0.401336339818801,0.401844497039815;
%     3,2.62405263065905,3.75813344832298,2.98371674540276,3.81799982863915,2.78747075469498,3.29471668335911];
% positions_for_A=[1,3.72102317542056,1.46803574249391,3.70674065758628,2.46242533721697,1.78982055662434,0.609623707440825;
%     0,1.47020128084852,0.612511298166949,0.812323004640796,0.545540103526116,0.199189444075264,0.397466995795309];
% 
no_A=length(positions_for_A(1,:));
no_B=length(positions_for_B(1,:));
no_C=length(positions_for_C(1,:));

figure(2)
for i=1:no_A
    while(positions_for_A(2,i) < 0) % Keep the objects above the surface
        positions_for_A(2,i)=normrnd(position_A(2),0.5,1,1);
    end
    while(positions_for_A(1,i) < 0) % Keep the objects above the surface
        positions_for_A(1,i)=normrnd(position_A(1),0.5,1,1);
    end
    while(positions_for_B(2,i) < 0)
        positions_for_B(2,i)=normrnd(position_B(2),0.5,1,1);
    end
    while(positions_for_B(1,i) < 0)
        positions_for_B(1,i)=normrnd(position_B(1),0.5,1,1);
    end
    while(positions_for_C(2,i) < 0)
        positions_for_C(2,i)=normrnd(position_C(2),0.5,1,1);
    end
    while(positions_for_C(1,i) < 0)
        positions_for_C(1,i)=normrnd(position_C(1),0.5,1,1);
    end
    rectangle('Position',[positions_for_A(1,i) positions_for_A(2,i) 2 3],'EdgeColor','b','LineWidth',3);
    hold on
    rectangle('Position',[positions_for_B(1,i) positions_for_B(2,i) 4 2],'EdgeColor','r','LineWidth',3);
    hold on
    rectangle('Position',[positions_for_C(1,i) positions_for_C(2,i) 2 2],'EdgeColor','g','LineWidth',3);
    hold on
    axis equal
    axis([-10 10 0 10]);
end
hold off

%Assign unary potentials from observtions--uniform unary potentials
unaryA=ones(no_A,1);
unaryB=ones(no_B,1);
unaryC=ones(no_C,1);

% for i=1:length(positions_for_A(1,:))
%     sigma1 = 0.1;
%     sigma2 = 0.5;
%     if(positions_for_A(2,i) > 1.5)
%         unaryA(:,i)= exp(-dist(positions_for_A(2,i))/2*sigma1^2);
%     elseif (positions_for_A(2,i) <1)
%         unaryA(:,i)=[0.5;0.5];
%     end
%     if(positions_for_B(2,i) > 1.5)
%         unaryB(:,i)=[0.4;0.6];
%     elseif (positions_for_B(2,i) <1)
%         unaryB(:,i)=[0.5;0.5];
%     end
% 
% end

%Assign Edge Potentials . Only a pair from (A,B) have edge potenitals not
%(A,A) and not (B,B)

BinaryAB= zeros(no_A,no_B);
Sigma=0.04;

for i=1:no_A
    for j=1:no_B
        %calculate the approximate center of the bounding boxes using the
        % posiitons and dimensions of the boxes
        x1= positions_for_A (1,i) +  w1/2; % I do not quite know if the objects will be in the directions of axes
        y1= positions_for_A (2,i) +  h1/2;
        x2= positions_for_B (1,j) +  w2/2; % I do not quite know if the objects will be in the directions of axes
        y2= positions_for_B (2,j) +  h2/2;
        if (x1>=x2)
            if (y1>=y2)
                if (y1==y2)
                    BinaryAB(i,j)=0.8;
                elseif((y1-y2)>1.2*((h1+h2)/2))
                    if((x1-x2)>1.1*(w1+w2)/2)
                      %Binary(i,j) = exp(-(((x1-x2)-((w1+w2)/2))^2+ ((y1-y2)-((h1+h2)/2))^2)/(2*Sigma));
                      BinaryAB(i,j) = 0.001;%exp(-(min((x1-x2),(y1-y2))^2)/(2*Sigma));
                    elseif((x1-x2)<0.9*(w1+w2)/2)
                      BinaryAB(i,j)= 0.01;   
                    end
                elseif((y1-y2)<0.8*((h1+h2)/2))
                    if((x1-x2)>= 1.1 *(w1+w2)/2)
                      BinaryAB(i,j) = 0.001%0.8* exp(-(y1-h1/2)^2-2*(y2-h2/2)^2/(2*Sigma));
                    elseif((x1-x2)< 0.9 *(w1+w2)/2)
                      BinaryAB(i,j)= 0.01;
                    end
                else
                    if((x1-x2)==(w1+w2)/2 || (y1-y2)==((h1+h2)/2))
                       BinaryAB(i,j)=0.8;
                    end
                end
            else
                if((y2-y1)>=1.2*((h1+h2)/2))
                    if((x1-x2)>=1.2*(w1+w2)/2)
                     % Binary(i,j) = exp(-(((x1-x2)-((w1+w2)/2))^2+ ((y2-y1)-((h1+h2)/2))^2)/(2*Sigma));
                      BinaryAB(i,j) = 0.01 ;%exp(-(min((x1-x2),(y2-y1))^2)/(2*Sigma));
                    elseif ((x1-x2)<0.9*(w1+w2)/2)
                      BinaryAB(i,j)= 0.01;% exp(-(x1-x2)^2/(2*Sigma));
                    end
                elseif((y2-y1)<0.8*((h1+h2)/2))
                    if((x1-x2)>=1.1 * (w1+w2)/2)
                      BinaryAB(i,j) =0.001; %0.8* exp(-(y2-h2/2)^2-2*(y1-h1/2)^2/(2*Sigma));
                    elseif((x1-x2)< 0.9 * (w1+w2)/2)
                      BinaryAB(i,j)= 0.01;
                    end
                else
                     if((x1-x2)==(w1+w2)/2 || (y2-y1)==((h1+h2)/2))
                         BinaryAB(i,j)=0.8;
                     end
                end
            end
        else
            if (y1>=y2)
                if (y1==y2)
                    BinaryAB(i,j)=0.8;
                elseif((y1-y2)>1.2 *((h1+h2)/2))
                    if((x2-x1)>=1.1*(w1+w2)/2)
                      %Binary(i,j) = exp(-(((x2-x1)-((w1+w2)/2))^2+ ((y1-y2)-((h1+h2)/2))^2)/(2*Sigma));
                      BinaryAB(i,j) = 0.01; %exp(-(min((x2-x1),(y1-y2))^2)/(2*Sigma));

                    elseif((x2-x1)<0.9*(w1+w2)/2)
                      BinaryAB(i,j)= 0.01; %exp(-(x2-x1)^2/(2*Sigma));
                    end
                elseif((y1-y2)<0.8 *((h1+h2)/2))
                    if((x2-x1)>=1.1*(w1+w2)/2)
                      BinaryAB(i,j) = 0.001; %0.8 * exp(-(y1-h1/2)^2-2*(y2-h2/2)^2/(2*Sigma));
                    elseif((x2-x1)<0.9*(w1+w2)/2)
                      BinaryAB(i,j)= 0.01;
                    end
                else
                      if ((x2-x1)==(w1+w2)/2 && (y1-y2)==((h1+h2)/2))
                          BinaryAB(i,j)=0.8;
                      end
                end
            else
                if((y2-y1)>=1.2*((h1+h2)/2))
                    if((x2-x1)>=1.2*(w1+w2)/2)
%                       Binary(i,j) = exp(-(((x2-x1)-((w1+w2)/2))^2+ ((y2-y1)-((h1+h2)/2))^2)/(2*Sigma));
                      BinaryAB(i,j) = 0.01 ;%exp(-(min((x2-x1),(y2-y1))^2)/(2*Sigma));

                    elseif((x2-x1)<0.9*(w1+w2)/2)
                      BinaryAB(i,j)= 0.01; %exp(-(x1-x2)^2/(2*Sigma));
                    end
                elseif((y2-y1)<0.9*((h1+h2)/2))
                    if((x2-x1)>=1.1*(w1+w2)/2)
                      BinaryAB(i,j) = 0.001; %0.8 * exp(-(y2-h2/2)^2-2*(y1-h1/2)^2/(2*Sigma));
                    elseif((x2-x1)<0.9*(w1+w2)/2)
                      BinaryAB(i,j)= 0.01;
                    end
                else
                    if(((x2-x1)==(w1+w2)/2) && (y2-y1)==((h1+h2)/2))
                       BinaryAB(i,j)=0.8;
                    end
                end
            end
        end
    end
end

BinaryAB=BinaryAB/sum(sum(BinaryAB));

BinaryAC= zeros(no_A,no_C);
Sigma=0.04;

for i=1:no_A
    for j=1:no_C
        %calculate the approximate center of the bounding boxes using the
        % posiitons and dimensions of the boxes
        x1= positions_for_A (1,i) +  w1/2; % I do not quite know if the objects will be in the directions of axes
        y1= positions_for_A (2,i) +  h1/2;
        x2= positions_for_C (1,j) +  w2/2; % I do not quite know if the objects will be in the directions of axes
        y2= positions_for_C (2,j) +  h2/2;
        if (x1>=x2)
            if (y1>=y2)
                if (y1==y2)
                    BinaryAC(i,j)=0.8;
                elseif((y1-y2)>1.2*((h1+h2)/2))
                    if((x1-x2)>1.1*(w1+w2)/2)
                      %Binary(i,j) = exp(-(((x1-x2)-((w1+w2)/2))^2+ ((y1-y2)-((h1+h2)/2))^2)/(2*Sigma));
                      BinaryAC(i,j) = exp(-(min((x1-x2),(y1-y2))^2)/(2*Sigma));
                    elseif((x1-x2)<0.9*(w1+w2)/2)
                      BinaryAC(i,j)= 0.01;   
                    end
                elseif((y1-y2)<0.8*((h1+h2)/2))
                    if((x1-x2)>= 1.1 *(w1+w2)/2)
                      BinaryAC(i,j) = 0.001; %0.8* exp(-(y1-h1/2)^2-2*(y2-h2/2)^2/(2*Sigma));
                    elseif((x1-x2)< 0.9 *(w1+w2)/2)
                      BinaryAC(i,j)= 0.01;
                    end
                else
                    if((x1-x2)==(w1+w2)/2 || (y1-y2)==((h1+h2)/2))
                       BinaryAC(i,j)=0.8;
                    end
                end
            else
                if((y2-y1)>=1.2*((h1+h2)/2))
                    if((x1-x2)>=1.2*(w1+w2)/2)
                     % Binary(i,j) = exp(-(((x1-x2)-((w1+w2)/2))^2+ ((y2-y1)-((h1+h2)/2))^2)/(2*Sigma));
                      BinaryAC(i,j) = 0.01 ;%exp(-(min((x1-x2),(y2-y1))^2)/(2*Sigma));
                    elseif ((x1-x2)<0.9*(w1+w2)/2)
                      BinaryAC(i,j)= 0.01;% exp(-(x1-x2)^2/(2*Sigma));
                    end
                elseif((y2-y1)<0.8*((h1+h2)/2))
                    if((x1-x2)>=1.1 * (w1+w2)/2)
                      BinaryAC(i,j) =0.001 %0.8* exp(-(y2-h2/2)^2-2*(y1-h1/2)^2/(2*Sigma));
                    elseif((x1-x2)< 0.9 * (w1+w2)/2)
                      BinaryAC(i,j)= 0.01;
                    end
                else
                     if((x1-x2)==(w1+w2)/2 || (y2-y1)==((h1+h2)/2))
                         BinaryAC(i,j)=0.8;
                     end
                end
            end
        else
            if (y1>=y2)
                if (y1==y2)
                    BinaryAC(i,j)=0.8;
                elseif((y1-y2)>1.2 *((h1+h2)/2))
                    if((x2-x1)>=1.1*(w1+w2)/2)
                      %Binary(i,j) = exp(-(((x2-x1)-((w1+w2)/2))^2+ ((y1-y2)-((h1+h2)/2))^2)/(2*Sigma));
                      BinaryAC(i,j) = 0.01; %exp(-(min((x2-x1),(y1-y2))^2)/(2*Sigma));

                    elseif((x2-x1)<0.9*(w1+w2)/2)
                      BinaryAC(i,j)= 0.01; %exp(-(x2-x1)^2/(2*Sigma));
                    end
                elseif((y1-y2)<0.8 *((h1+h2)/2))
                    if((x2-x1)>=1.1*(w1+w2)/2)
                      BinaryAC(i,j) = 0.001 ;%0.8 * exp(-(y1-h1/2)^2-2*(y2-h2/2)^2/(2*Sigma));
                    elseif((x2-x1)<0.9*(w1+w2)/2)
                      BinaryAC(i,j)= 0.01;
                    end
                else
                      if ((x2-x1)==(w1+w2)/2 && (y1-y2)==((h1+h2)/2))
                          BinaryAC(i,j)=0.8;
                      end
                end
            else
                if((y2-y1)>=1.2*((h1+h2)/2))
                    if((x2-x1)>=1.2*(w1+w2)/2)
%                       Binary(i,j) = exp(-(((x2-x1)-((w1+w2)/2))^2+ ((y2-y1)-((h1+h2)/2))^2)/(2*Sigma));
                      BinaryAC(i,j) = 0.01 ;%exp(-(min((x2-x1),(y2-y1))^2)/(2*Sigma));

                    elseif((x2-x1)<0.9*(w1+w2)/2)
                      BinaryAC(i,j)= 0.01; %exp(-(x1-x2)^2/(2*Sigma));
                    end
                elseif((y2-y1)<0.9*((h1+h2)/2))
                    if((x2-x1)>=1.1*(w1+w2)/2)
                      BinaryAC(i,j) = 0.001 ;%0.8 * exp(-(y2-h2/2)^2-2*(y1-h1/2)^2/(2*Sigma));
                    elseif((x2-x1)<0.9*(w1+w2)/2)
                      BinaryAC(i,j)= 0.01;
                    end
                else
                    if(((x2-x1)==(w1+w2)/2) && (y2-y1)==((h1+h2)/2))
                       BinaryAC(i,j)=0.8;
                    end
                end
            end
        end
    end
end

BinaryAC=BinaryAC/sum(sum(BinaryAC));


BinaryCB= zeros(no_C,no_B);
Sigma=0.04;

for i=1:no_C
    for j=1:no_B
        %calculate the approximate center of the bounding boxes using the
        % posiitons and dimensions of the boxes
        x1= positions_for_C (1,i) +  w1/2; % I do not quite know if the objects will be in the directions of axes
        y1= positions_for_C (2,i) +  h1/2;
        x2= positions_for_B (1,j) +  w2/2; % I do not quite know if the objects will be in the directions of axes
        y2= positions_for_B (2,j) +  h2/2;
        if (x1>=x2)
            if (y1>=y2)
                if (y1==y2)
                    BinaryCB(i,j)=0.8;
                elseif((y1-y2)>1.2*((h1+h2)/2))
                    if((x1-x2)>1.1*(w1+w2)/2)
                      %Binary(i,j) = exp(-(((x1-x2)-((w1+w2)/2))^2+ ((y1-y2)-((h1+h2)/2))^2)/(2*Sigma));
                      BinaryCB(i,j) = 0.001; %exp(-(min((x1-x2),(y1-y2))^2)/(2*Sigma));
                    elseif((x1-x2)<0.9*(w1+w2)/2)
                      BinaryCB(i,j)= 0.01;   
                    end
                elseif((y1-y2)<0.8*((h1+h2)/2))
                    if((x1-x2)>= 1.1 *(w1+w2)/2)
                      BinaryCB(i,j) = 0.001 %0.8* exp(-(y1-h1/2)^2-2*(y2-h2/2)^2/(2*Sigma));
                    elseif((x1-x2)< 0.9 *(w1+w2)/2)
                      BinaryCB(i,j)= 0.01;
                    end
                else
                    if((x1-x2)==(w1+w2)/2 || (y1-y2)==((h1+h2)/2))
                       BinaryCB(i,j)=0.8;
                    end
                end
            else
                if((y2-y1)>=1.2*((h1+h2)/2))
                    if((x1-x2)>=1.2*(w1+w2)/2)
                     % Binary(i,j) = exp(-(((x1-x2)-((w1+w2)/2))^2+ ((y2-y1)-((h1+h2)/2))^2)/(2*Sigma));
                      BinaryCB(i,j) = 0.01 ;%exp(-(min((x1-x2),(y2-y1))^2)/(2*Sigma));
                    elseif ((x1-x2)<0.9*(w1+w2)/2)
                      BinaryCB(i,j)= 0.01;% exp(-(x1-x2)^2/(2*Sigma));
                    end
                elseif((y2-y1)<0.8*((h1+h2)/2))
                    if((x1-x2)>=1.1 * (w1+w2)/2)
                      BinaryCB(i,j) =0.001; %0.8* exp(-(y2-h2/2)^2-2*(y1-h1/2)^2/(2*Sigma));
                    elseif((x1-x2)< 0.9 * (w1+w2)/2)
                      BinaryCB(i,j)= 0.01;
                    end
                else
                     if((x1-x2)==(w1+w2)/2 || (y2-y1)==((h1+h2)/2))
                         BinaryCB(i,j)=0.8;
                     end
                end
            end
        else
            if (y1>=y2)
                if (y1==y2)
                    BinaryCB(i,j)=0.8;
                elseif((y1-y2)>1.2 *((h1+h2)/2))
                    if((x2-x1)>=1.1*(w1+w2)/2)
                      %Binary(i,j) = exp(-(((x2-x1)-((w1+w2)/2))^2+ ((y1-y2)-((h1+h2)/2))^2)/(2*Sigma));
                      BinaryCB(i,j) = 0.01; %exp(-(min((x2-x1),(y1-y2))^2)/(2*Sigma));

                    elseif((x2-x1)<0.9*(w1+w2)/2)
                      BinaryCB(i,j)= 0.01; %exp(-(x2-x1)^2/(2*Sigma));
                    end
                elseif((y1-y2)<0.8 *((h1+h2)/2))
                    if((x2-x1)>=1.1*(w1+w2)/2)
                      BinaryCB(i,j) = 0.001; %0.8 * exp(-(y1-h1/2)^2-2*(y2-h2/2)^2/(2*Sigma));
                    elseif((x2-x1)<0.9*(w1+w2)/2)
                      BinaryCB(i,j)= 0.01;
                    end
                else
                      if ((x2-x1)==(w1+w2)/2 && (y1-y2)==((h1+h2)/2))
                          BinaryCB(i,j)=0.8;
                      end
                end
            else
                if((y2-y1)>=1.2*((h1+h2)/2))
                    if((x2-x1)>=1.2*(w1+w2)/2)
%                       Binary(i,j) = exp(-(((x2-x1)-((w1+w2)/2))^2+ ((y2-y1)-((h1+h2)/2))^2)/(2*Sigma));
                      BinaryCB(i,j) = 0.01 ;%exp(-(min((x2-x1),(y2-y1))^2)/(2*Sigma));

                    elseif((x2-x1)<0.9*(w1+w2)/2)
                      BinaryCB(i,j)= 0.01; %exp(-(x1-x2)^2/(2*Sigma));
                    end
                elseif((y2-y1)<0.9*((h1+h2)/2))
                    if((x2-x1)>=1.1*(w1+w2)/2)
                      BinaryCB(i,j) = 0.001 ; %0.8 * exp(-(y2-h2/2)^2-2*(y1-h1/2)^2/(2*Sigma));
                    elseif((x2-x1)<0.9*(w1+w2)/2)
                      BinaryCB(i,j)= 0.01;
                    end
                else
                    if(((x2-x1)==(w1+w2)/2) && (y2-y1)==((h1+h2)/2))
                       BinaryCB(i,j)=0.8;
                    end
                end
            end
        end
    end
end

BinaryCB=BinaryCB/sum(sum(BinaryCB));


%Normalize the Binary potentials
% sum=0;
%  for i=1:no_A
%      sum=sum+sum_a(Binary(i,:));
%  end
%  
%  Binary=Binary/sum;
  

%Initialize previous messages

prev_messagefromAtoB = ones(no_A,1);
prev_messagefromAtoC = ones(no_A,1);
prev_messagefromBtoA = ones(no_B, 1);
prev_messagefromBtoC = ones(no_B, 1);
prev_messagefromCtoA = ones(no_C, 1);
prev_messagefromCtoB = ones(no_C, 1);


%Initialize the beliefs

belA=zeros(no_A,1);
belB=zeros(no_B,1);
belC=zeros(no_C,1);

%Now start the loop- There is only one edge

for iter=15
    
    % Calculate the messages to all B nodes from all A nodes , I don't think
    % messages are commutative in this case
    
         messagefromAtoB= BinaryAB*unaryA.*prev_messagefromCtoA;
         messagefromBtoA= BinaryAB'*unaryB.*prev_messagefromCtoB;
         messagefromAtoC= BinaryAC*unaryA.*prev_messagefromBtoA;
         messagefromBtoC= BinaryCB'*unaryB.*prev_messagefromAtoB;
         messagefromCtoA= BinaryAC'*unaryC.*prev_messagefromBtoC;
         messagefromCtoB= BinaryCB*unaryC.*prev_messagefromAtoC;
         
         belA=unaryA.*messagefromBtoA.*messagefromCtoA;
         belB=unaryB.*messagefromAtoB.*messagefromCtoB;
         belC=unaryC.*messagefromAtoC.*messagefromBtoC;
         
         %Normalize the beliefs
         belA=belA/sum_a(belA);
         belB=belB/sum_a(belB);
         belC=belB/sum_a(belC);
         
         %update the messages
         prev_messagefromBtoA=messagefromBtoA;
         prev_messagefromAtoB=messagefromAtoB;
         prev_messagefromCtoA=messagefromCtoA;
         prev_messagefromBtoC=messagefromBtoC;
         prev_messagefromAtoC=messagefromAtoC;
         prev_messagefromCtoB=messagefromCtoB;
   
end
max_A=max(belA);
max_B=max(belB);
max_C=max(belC);


for i=1:no_A
    if (max_A==belA(i))
        index_A=i;
        break;
    end
end

for i=1:no_B
    if (max_B==belB(i))
        index_B=i;
        break;
    end
end
for i=1:no_C
    if (max_C==belC(i))
        index_C=i;
        break;
    end
end

        
figure(3)
rectangle('Position',[positions_for_A(1,index_A) positions_for_A(2,index_A) w1 h1],'EdgeColor','b','LineWidth',3);
    hold on
rectangle('Position',[positions_for_B(1,index_B) positions_for_B(2,index_B) w2 h2],'EdgeColor','r','LineWidth',3);
    hold on
rectangle('Position',[positions_for_C(1,index_C) positions_for_C(2,index_C) w3 h3],'EdgeColor','g','LineWidth',3);
    hold on
    axis equal
    axis([-10 10 0 10]);












