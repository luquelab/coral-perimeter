function d = richardsonDistance(X,Y,len,step_max)
% initialize as a complex vector
z=X+1i*Y;
% count the number of elements
N=length(z);
% initialize the reference index and distance travelled
i0=1;
d=0;
while(i0<N)
   % create a temporary index vector to ensure we only look forward
   pmax=min(N,i0+step_max);
   idx_tmp=(i0+1):pmax;
   % calculate the euclidean distance of each point relative to the current
   % index
   ztmp=z(idx_tmp)-z(i0);
   D=abs(ztmp);
   % find all points within the euclidean distance "len"
   logik=(D>len);
   % this statement tells us to exit the loop whenever there are no points
   % beyond the radius, and add the "partial length" to the number of 
   % lengths travelled
   if(sum(logik)==0) % if there are no zeros in the logik vector
       i0=N; % put us at the end so we exit the while
       d=d+D(end)/len; % add the partial distance to the output var
   else
       % otherwise, find the first zero; the previous point is the boundary
       ii=find(logik,1,'first');
       % set the new index
       i0=idx_tmp(ii);
       % add a length to the output distance "d"
       d=d+1;
   end
end

% d=d/len;

end