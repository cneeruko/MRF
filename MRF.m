nColumns = 3;
nRows = 3;
nNodes = nColumns*nRows;
nStates = 2;
adj = zeros(nNodes);
bel=zeros(nNodes,2);
% 
% bel=ones(nNodes,1);
for r = 1:nRows
	for s = 1:nColumns
		if s < nColumns
			adj(s + nColumns*(r-1),(s+1) + nColumns*(r-1)) = 1; % Passenger behind
		end
		if r < nRows
			adj(s + nColumns*(r-1),s + nColumns*(r)) = 1; % Passenger to the right
		end
	end
end
adj = adj+adj';
edgeStruct = UGM_makeEdgeStruct(adj,nStates,0);
edgeStruct.edgeEnds=double(edgeStruct.edgeEnds);
%Define Node potenitals and Edge Potentials
for i=1:nNodes
    if((i== 1)||(i== 5) ||(i== 9))
        nodePot(:,i)=[0.3 ;
                        0.7 ];
    else
        nodePot(:,i)=[0.6 ;
                        0.4 ];
    end
end
for i=1:edgeStruct.nEdges
    if (i== 1)||(i==2)||(i==4)||(i==6)||(i==8)||(i==9)||(i==10)
       edgePot(:,:,i)=[0.3 0.7;
                       0.7 0.3] ;
    else
        edgePot(:,:,i)=[0.6 0.4;
                       0.4 0.6] ;
    end
 
end
param=ones(2*edgeStruct.nEdges,2);
edgeStruct.edgeEnds(edgeStruct.nEdges+1:2*edgeStruct.nEdges,1)=edgeStruct.edgeEnds(1:edgeStruct.nEdges,2);
edgeStruct.edgeEnds(edgeStruct.nEdges+1:2*edgeStruct.nEdges,2)=edgeStruct.edgeEnds(1:edgeStruct.nEdges,1);
edgeStruct.edgeEnds=[edgeStruct.edgeEnds param];

x= [1:4];
y = [1:4];
[X,Y] = meshgrid(x,y);
Z= [reshape(bel(:,2),3,3) zeros(3,1);
    zeros(1,4)];
set(1,'renderer','zbuffer');  % appears to work better for recording movie
    vidObj = VideoWriter('APsim','MPEG-4');
     set(vidObj, 'FrameRate', 3); %min(5, 1/pauseLen)
    open(vidObj);
figure(1)
colormap('parula(10)');
h=pcolor(X,Y,Z);
set(gca,'clim',[0 1]);
currFrame = getframe(gca);
writeVideo(vidObj,currFrame);

%loop
for iter=1:20 %edgeStruct.maxIter
       
  for edge=1:2*edgeStruct.nEdges
    s= edgeStruct.edgeEnds(edge,1);
    t= edgeStruct.edgeEnds(edge,2);
   
    %for each t, product of messages from all neighbors to s except t
%     tnbrs=find(adj(t,:));
    p=[1;1];
    row=0;
    sum=0;
%     for i=1:length(tnbrs) %all s nodes
%         s=tnbrs(i);
        %product term
        snbrs=find(adj(s,:));
        for j=1:length(snbrs)
            u=snbrs(j);
            if(t~=u)
              row=edge_row(u,s,edgeStruct.edgeEnds);
              if(row==0)
                  exit;
              else 
                 msg=double(edgeStruct.edgeEnds(row,1+2*iter:2+2*iter))';
                 p=p.*msg;
              end
            end
        end 
        strow=edge_row(s,t,edgeStruct.edgeEnds);
        if(rem(strow,12)==0)
         p=edgePot(:,:,strow/12)*(p.*nodePot(:,s));
        else
         p=edgePot(:,:,rem(strow,12))*(p.*nodePot(:,s));
        end
         edgeStruct.edgeEnds(strow,3+2*iter:4+2*iter)=p;
%          sum=sum+p;
       
%     end
%     node_mes(t)= sum;
     
  end
  
   %For each s, Caculate its belief in each iteration
  for s=1:nNodes
      p=[1;1];
      snbrs=find(adj(s,:));
      for j=1:length(snbrs)
          t=snbrs(j);
          row=edge_row(t,s,edgeStruct.edgeEnds);
          if(row==0)
             exit;
          else
             p=p.*edgeStruct.edgeEnds(row,3+2*iter:4+2*iter)';
          end
      end
      if (p(1)==0 && p(2)==0)
          bel(s,:)=round(bel(s,:));
      else
      bel(s,:)=(nodePot(:,s).*p)';
      end
      %Normalize belief
      bel(s,:)=bel(s,:)/sum_a(bel(s,:)');
      
      product(s)= bel(s,1)*bel(s,2);
  end
  %visualize
  
  Z= [reshape(bel(:,2),3,3) zeros(3,1);
    zeros(1,4)];
set(h,'xdata',X,'ydata',Y,'cdata',Z);
axis ij
title(['Iteration :',num2str(iter)]);
currFrame = getframe(gca);
writeVideo(vidObj,currFrame);
pause(1);

if (length(find(product))==0)
    title('Model converged');
    break;
end

po=0;
  
end
close(vidObj);
        vidObj = [];

