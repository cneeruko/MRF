function  output = mrf_flower( path )
% karthik: given a image this asks for foreground and then does bg subtraction 
% based on mrf belief propagation

% yi be the RGB values of the image
% xi be the classification {0, 1} i.e. 1 meaning its the foreground
% We need to compute the joint distribution 
% p({x}, {y}) = 1/z prod_ij(psi_ij(xi, xj)) prod_i(phi_i(xi, yi))
% psi_ij -> compatibility function will call pairwise(i, j)
% phi_i -> joint compatibility function unary(i)

% Define functions pairwise(i, j)
% For pairwise each of the pixel has 4 neighbours so lets have RXCX4
% matrix
% Define function unary(i)
% For unary term there is value for each pixel so lets have a mat RXC

% First load the image
[I,map] = imread(path);
I = double(I);
r = size(I, 1);
c = size(I, 2);

L = zeros(r, c);

 % Set a fg rgb value
fg = [230;173;84];
bg = [8;31;3];
 
%Initialize the belief of all the pixels
belM_f = zeros(r, c,2);
%belM_b = zeros(r, c);
%Each pixel will have a 2-by-1 node potenntial which will be derived from
%observation
unaryM = 0.5*ones(r, c,2); 
%Each edge potential will be a 2*2 matrix 
pairM = 0.5*ones(r, c, 4, 2,2);  

 prev_messageM = ones(r, c, 5, 2);
 messageM = ones(r, c, 5, 2);
% msgM = ones(r, c, 5);

 change = ones(r, c);
 
 % update unary Potential based on the RGB values
 [unaryM, L] = updateUnary(unaryM, I, fg, bg, L);
%  x= [1:300];
%  y = [1:199];
% [X,Y] = meshgrid(x,y);
  figure(1)
%  map=colormap(map);
%  h=pcolor(X,Y,L);
%  set(h, 'EdgeColor', 'none');
%  axis ij
 h=imshow(L,map);
 set(gca,'clim',[0 1]);
 

 pairM = updateBinary(pairM, unaryM, L);
 
 %Loop untill belief doesn't change
 count =1;
 
 while(size(size(change~=0))~=0)
   prev_B = belM_f;
    for i=2:r-1
        for j=2:c-1
             %Send message on each edge of a pixel
             %update the messageM(i, j)
             
            %For foreground
             % For target i, j message coming from left->node
            messageM(i, j, 2, :) = reshape(pairM(i, j-1, 1, :,:),2,2)*(reshape(unaryM(i, j-1, :),2,1).*reshape(prev_messageM(i, j-1, 1, :),2,1).*reshape(prev_messageM(i, j-1, 3, :),2,1).*reshape(prev_messageM(i, j-1, 4, :),2,1));
             % For target i, j message coming from right->node
            messageM(i, j, 1,:) = reshape(pairM(i, j+1, 2, :,:),2,2)*(reshape(unaryM(i, j+1, :),2,1).*reshape(prev_messageM(i, j+1, 2, :),2,1).*reshape(prev_messageM(i, j+1, 3, :),2,1).*reshape(prev_messageM(i, j+1, 4, :),2,1));
            % For target i, j message coming from up->node
            messageM(i, j, 4, :) = reshape(pairM(i-1, j, 3, :,:),2,2)*(reshape(unaryM(i-1, j, :),2,1).*reshape(prev_messageM(i-1, j, 1, :),2,1).*reshape(prev_messageM(i-1, j, 2, :),2,1).*reshape(prev_messageM(i-1, j, 3, :),2,1));
            % For target i, j message coming from down->node
            messageM(i, j, 3, :) = reshape(pairM(i+1, j, 4, :,:),2,2)*(reshape(unaryM(i+1, j, :),2,1).*reshape(prev_messageM(i+1, j, 1, :),2,1).*reshape(prev_messageM(i+1, j, 2, :),2,1).*reshape(prev_messageM(i+1, j, 4, :),2,1));
            %messageM(i, j, 5, 1) = messageM(i, j, 1, 1)+messageM(i, j, 2, 1)+messageM(i, j, 3, 1)+messageM(i, j, 4, 1);

%             %For background
%             % For target i, j message coming from left->node
%             messageM(i, j, 1, 2) = unaryM(i, j-1, 2)*pairM(i, j-1, 1, 2)*prev_messageM(i, j-1, 1, 2)*prev_messageM(i, j-1, 3, 2)*prev_messageM(i, j-1, 4, 2);
%             % For target i, j message coming from right->node
%             messageM(i, j, 2, 2) = unaryM(i, j+1, 2)*pairM(i, j+1, 2, 2)*prev_messageM(i, j+1, 2, 2)*prev_messageM(i, j+1, 3, 2)*prev_messageM(i, j+1, 4, 2);
%             % For target i, j message coming from up->node
%             messageM(i, j, 3, 2) = unaryM(i-1, j, 2)*pairM(i-1, j, 3, 2)*prev_messageM(i-1, j, 1, 2)*prev_messageM(i-1, j, 2, 2)*prev_messageM(i-1, j, 3, 2);
%             % For target i, j message coming from down->node
%             messageM(i, j, 4, 2) = unaryM(i+1, j, 2)*pairM(i+1, j, 4, 2)*prev_messageM(i+1, j, 1, 2)*prev_messageM(i+1, j, 2, 2)*prev_messageM(i+1, j, 4, 2);
%             %messageM(i, j, 5, 2) = messageM(i, j, 1, 2)+messageM(i, j, 2, 2)+messageM(i, j, 3, 2)+messageM(i, j, 4, 2);
             
%             for a=1:4
%                 msgM(i, j, a) = messageM(i, j, a, 1) + messageM(i, j, a, 2);
%             end  
        end
    end
    k=1;
   for i=1:r
       for j=1:c
           
         if (belM_f(i,j,1) > 0.7)  
               continue;
          else
           %update the belM(i, j)
           %belM_f(i, j) = unaryM(i, j, 1)*messageM(i, j, 1, 1)*messageM(i, j, 2, 1)*messageM(i, j, 3, 1)*messageM(i, j, 4, 1);
           %belM_b(i, j) = unaryM(i, j, 2)*messageM(i, j, 1, 2)*messageM(i, j, 2, 2)*messageM(i, j, 3, 2)*messageM(i, j, 4, 2); 
           belM_f(i, j,:) = reshape(unaryM(i, j, :),2,1).*reshape(messageM(i, j, 1, :),2,1).*reshape(messageM(i, j, 2, :),2,1).*reshape(messageM(i, j, 3, :),2,1).*reshape(messageM(i, j, 4, :),2,1);
          % belM_b(i, j) = unaryM(i, j, 2)*msgM(i, j, 1)*msgM(i, j, 2)*msgM(i, j, 3)*msgM(i, j, 4);
           
           Nor = belM_f(i, j,1) + belM_f(i, j,2);
           belM_f(i, j,:) = belM_f(i, j,:)/Nor;
          % belM_b(i, j) = belM_b(i, j,2)/Nor;
          
           product(k)= belM_f(i,j,1)*belM_f(i,j,2);
           k=k+1;
         end
       end
   end
    prev_messageM = messageM;
    next_B = belM_f;
    change = abs(next_B - prev_B);    
    temp_image = belM_f(:,:,1); 
    set(h,'cdata',temp_image);
%     axis ij
    if (length(find(product))==0)
      title('Model converged');
      break;
    end
    title(['Iteration :',num2str(count)]);
    %imshow(temp_image);
    count = count+1
    pause(0.5);
    
 end
 
 end
 
function [unaryM, L] = updateUnary(M, I, fg, bg, L)
    % Assume a variance for the prob generation
    sigma1 = 0.1;
    sigma2 = 0.5;
    for i=1:size(M, 1)
        for j=1:size(M, 2)
            M(i, j,1) = exp(-dist(I(i, j, :), fg)/2*sigma1^2)  ;
            M(i, j,2) = exp(-dist(I(i, j, :), bg)/2*sigma2^2) ;
            if(M(i, j,1)>0.35)
%                 disp('yes');
                L(i, j)=1;
%                 M(i, j,1)=M(i, j,1)+0.7/(M(i, j,1)+M(i, j,2))+1;
%                 M(i, j,2)=M(i, j,1)+0.3/(M(i, j,1)+M(i, j,2))+1;
%             else
            %Normalize unary Potentials
%                 M(i, j,1)=M(i, j,1)+0.3/(M(i, j,1)+M(i, j,2))+1;
%                 M(i, j,2)=M(i, j,1)+0.7/(M(i, j,1)+M(i, j,2))+1;
            end
                M(i, j,:)=M(i, j,:)/(M(i, j,1)+M(i, j,2));

        end
    end
    unaryM = M;
    
end

function pairM = updateBinary(M, U, L)
    %Giving weights because of the potentials
    ws = [0.7 0.3;0.3 0.7];
    wd = [0.3 0.7;0.7 0.3];

    for i=2:size(U, 1)-1
        for j=2:size(U, 2)-1
            %left
            if(L(i, j)==L(i, j-1)) %Same label
                M(i, j, 1, :,:) = ws;%U(i, j, 1)*ws;
            else
                M(i, j, 1, :,:) = wd;%U(i, j, 1)*wd;
            end
            %right
            if(L(i, j)==L(i, j+1)) %Same label
                M(i, j, 2, :,:) = ws;%U(i, j, 1)*ws;
            else
                M(i, j, 2, :,:) = wd;%U(i, j, 1)*wd;
            end
            %up
            if(L(i, j)==L(i-1, j)) %Same label
                M(i, j, 3, :,:) = ws;%U(i, j, 1)*ws;
            else
                M(i, j, 3, :,:) = wd;%U(i, j, 1)*wd;
            end
            %down
            if(L(i, j)==L(i+1, j)) %Same label
                M(i, j, 4, :,:) = ws;%U(i, j, 1)*ws;
            else
                M(i, j, 4, :,:) = wd;%U(i, j, 1)*wd;
            end    
        end
    end
    pairM = M;

 end


function d = dist(A, B)
     d = ((A(:, :, 1)-B(1, 1))^2 + (A(:, :, 2)-B(2, 1))^2 +(A(:, :, 3)-B(3, 1))^2);
 end

