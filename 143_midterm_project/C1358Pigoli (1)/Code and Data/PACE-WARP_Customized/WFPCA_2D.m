% ===========
% Description:
% ===========
%            This program is a flavour of the original WFPCA.m in order to
%            use full matrix data rather than just arrays. It has limited
%            functionality compared to the original; namely it can not work
%            on sparse data. The rest of this description shows what the
%            original WPFA did and the current WFPCA_2D can not do.
%
%            Originally this program was for aligning warped functional
%            data. It used  a pairwise warping method to obtain the desired
%            alignment (registration) of the random trajectories. The data
%            could be regular (recorded on a regular grid) or sparse. PACE 
%            would be used as a preliminary step to evaluate the data on a 
%            regular time grid but here in the 2D case you HAVE TO HAVE A 
%            REGULAR GRID BEFORE YOU START THE ANALYSIS. The original
%            algorithm gave the aligned curves and the warping function; 
%            this algorithm gives out the aligned surfaces and the warping 
%            functions. Users could originally also choose to run PACE on 
%            the aligned curves but in this version they can't.
% 
%            We extend the idea that the model is Y(t)=X(h^{-1}(t))+error 
%            where Y(t) is the observed noisy curve/surface, X(t) is the 
%            aligned curve/surface and h^{-1}(t) is the inverse warping
%            function. The warping function h(t) can be used to align objects,
%            i.e. Y*(t)=Y(h(t)) is output aligned (still with noise, as 
%            only the measurement times are shifted according to the warping).
% 
%            References: Tang, R., M\"{u}ller, H.G. (2008). Pairwise curve
%            synchronization for functional data. Biometrika 95, 875-889.
%            Tang, R., M\"{u}ller, H.G. (2009). Time-synchronized 
%            clustering of gene expression trajectories. Biostatistics 10, 32-45.
%  
% ======
% Usage:
% ======
% 
% function [X,aligned,h,hinv,t_reg] = WFPCA_2D(y,t,nknots,lambda,choice)
%            Notice that unlike WPFA you can't have PACE output also.
% ======
% Input: 
% ======
%      y:          cell array of u n*m matrices; make sure that is on regular
%                  grid.
%      t:          1*m vector (if y is a matrix)
%      nknots:     number of knots used for estimating pairwise warping functions; 
%                  default value is 3.
%      lambda:     penalty parameter used for estimating pairwise warping
%                  functions (larger lambda will force the warping functions
%                  closer to the identity and thus result in less warping);
%                  use 0 as a first default value, which will set the penalty to V*10^-4, 
%                  where V is the average L2 norm of y-mean(y); otherwise
%                  provide lambda yourself: no fancy grid searches here as
%                  in the original WPCA.
%      choice:     method for estimating warping functions.         
%                  'weighted': use weighted averages of pairwise warping functions and choose 
%                              weights based on inverse pairwise distances  [default].         
%                  'truncated': pairs with large distances (the largest 10%) are truncated and the 
%                               averages of the remaining pairwise distances are used.         
% 
% =======
% Output:  
% =======  
%      X:         A cell array containing all returned values from PCA.m It is
%                 is always evaluated to be []. I keep it for the sake of familiarity. 
%      aligned:   cell aray with u n*m matrices where aligned{i} is the estimated aligned 
%                 surface evaluated at time vector t (for regular data)  for the ith subject.  
%      h:         n*m matrix where h(i,:) is the estimated warping function evaluated 
%                 at time vector t (for regular data) or t_reg (for sparse data) 
%                 for the ith subject (see below or Tang and M\"uller, 2008 for the definition).
%      hinv:      n*m matrix where hinv(i,:) is the estimated inverse warping function evaluated at 
%                 time vector t (for regular data) or t_reg (for sparse data) for the ith subject.
%      t_reg:     1*m vector, time grid for y_reg.
% 
%  
% 
%   To get individual values from X, type getVal(X,varname).
%   To obtain the names for X, type names(X).
%   To see an example, check example_WFPCA.m


function [X,aligned,h,hinv,t_reg] = WFPCA_2D(y,t,nknots,lambda,choice)
% A couple of arguments are missing compared to the usual WFPCA check above
% for details explaining what is missing

if nargin<5|isempty(choice) choice='weighted'; end
if nargin<4|isempty(lambda) lambda=0; end            
if nargin<3|isempty(nknots) nknots=3; end

u=length(y);
m=length(t);
[n,m] = size(y{1});

y_reg=y;
t_reg=t;

y_normalized={};

% supernorm normalization on 2-D
for i=1:u
    Z = y_reg{i};
    tempsca=max(abs(Z'));
    Z_normalized = Z.*repmat(1./tempsca',[1 m]); 
    y_normalized{i} = Z_normalized;
end 
clear tempsca;

yfd={};
% transfer normalized data into spline form
for i=1:u
     Z = y_normalized{i};
     for j=1:n
        Z_yfd(j,:)=spapi(2,t_reg,Z(j,:));  % transfering normalized data into spline form
     end 
     yfd{i} = Z_yfd;    
end

yfd_old={};
% transfer normalized data into spline form
for i=1:u
     Z = y_reg{i};
     for j=1:n
        Z_yfd_old(j)=spapi(2,t_reg,Z(j,:));  % transfering normalized data into spline form
     end 
     yfd_old{i} = Z_yfd_old;    
end

% use subset to estimate warping functions
subN=floor(min(max(30,sqrt(u)),u));  % Changed maximum subset size from 20 to 30.!!

% choose lambda if not specified
if ( isempty(lambda) || lambda==0 )
    %Concatenate all the 2-D variables in a huge matrix
    Z = zeros(u*n,m);
    for i=1:u
        Z(( (-1+i)*n+1:i*n),1:m) =  y_normalized{i};
    end
    Vy=sqrt(sum(trapz(t_reg',(Z-repmat(mean(Z),[(u*n),1]))'.^2))/((u*n)-1));
    lambda=Vy*10^-4;
    fprintf(1,['Default lambda value ' num2str(lambda) ' is used. \n']);
end        

indd=1;
hik=[]; 

for k=1:u
     tempind=randperm(u);
     tempind=tempind(1:subN);
     for i=1:subN;
          objecti=yfd{tempind(i)}; %Query curve
          objectk=yfd{k};          %Reference Curve
          if tempind(i)~=k
              %Get splines to matrix form
              ObjectImatrix = zeros(n,m);
              OjbectKmatrix = zeros(n,m);
              for p=1:1:n
                  ObjectImatrix(p,:) = objecti(p).coefs;
                  ObjectKmatrix(p,:) = objectk(p).coefs;
              end                 
              %%%Get the universal pairwise distance for all matrix rows
              
              %Get the mapping hik and the cost associated with it
              [hik(indd,:),ftemp]=rthik_E_2D( ObjectImatrix ,ObjectKmatrix ,t_reg,nknots,lambda);  
              %Get the warped object in a matrix form
              warped_object = zeros(n,m);
              for p=1:1:n
                   warped_object(p,:) = rtYH(t_reg,objecti(p),hik(indd,:));
              end  
              %Find the pairwise distance between the two curves
              %afterdistance(k,i)=mean((rtYH(t_reg,objecti,hik(indd,:))-fnval(objectk,t_reg)).^2);                                 
              afterdistance(k,i) = mean2((warped_object  - ObjectKmatrix ).^2);
              %"Jargon" for the averaging step
              lochat(indd,:)=[k,tempind(i)]; 
              indd=indd+1; 
          else  %so in case we are trying to warp a function with itself
              hik(indd,:)=t_reg;
              afterdistance(k,i)=0; 
              lochat(indd,:)=[k,k];  
              indd=indd+1;
          end    
     end
     clear tempind;
end

% take truncated averages or weighted averages of the pairwise warping functions 
h=[];
hinv=[];
aligned=[];
for k=1:u
    weight=[];
    if strcmp(choice,'weighted')
         tempw=afterdistance(k,:); 
         for i=1:subN
             if tempw(i)==0;
                 tempw(i)=inf;
             end
         end
         weight=1./tempw;
         clear tempw;
    elseif strcmp(choice,'truncated')
         weight=afterdistance(k,:)<=quantile(afterdistance(afterdistance>0),.90);            
    end
    if sum(weight)==0
        hinv(k,:)=t_reg;
    else
        weight=weight/sum(weight);
        hinv(k,:)=weight*hik(lochat(:,1)==k,:);
        hinv(k,m)=max(t_reg);
        hinv(k,1)=min(t_reg);
    end
    h(k,:)=fnval(spapi(2,hinv(k,:),t_reg),t_reg);
    h(k,1)=min(t_reg);
    h(k,m)=max(t_reg);
    %aligned(k,:)=fnval(yfd_old(k,:),h(k,:));  %Originally
    aligned_object= zeros(n,m);
    for i=1:1:p
        aligned_object(i,:) = fnval( yfd_old{k}(i) ,h(k,:) );
    end
    aligned{k} = aligned_object; %Save the aligned matrix in the cell array
end

X=[];

end
    
