function index = edge_row(u,s,edgeEnds)
 %compare the start and end nodes of the edge and give the corresponding
 %row number
 for i=1:size(edgeEnds,1)
     if(edgeEnds(i,1)==u && edgeEnds(i,2)==s)
         index=i;
         return
     else 
         index=0;
     end
 end
end