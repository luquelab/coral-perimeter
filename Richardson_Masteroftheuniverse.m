clear;
% read pix into str
% pix=dir('*.jpg');
pix=importdata('pc_values.csv',',',0);
%
for idx=1:length(pix.data)
    tstart=tic;
    if(strcmp(pix.textdata{idx},'CUR9_stitched_outlined_calibrated_220_BW.jpg'))
        continue;
    end
    display(pix.textdata{idx});
    I=imread(pix.textdata{idx});%'CUR8_alignment_cropped_caibrated_181_RGB_BW.jpg');
    pc = pix.data(idx);%181; %pixels/cm (resolution)
    BW=im2bw(I,.07);
    [rows,cols]=size(BW);
    %---------------------------------
    %find start and finish border pts along image edge
start=[];  fin=[];% define empty values for start and stop point

Left=BW(:,1);
Right=BW(:,cols);
Top=BW(1,:);
Bottom=BW(rows,:);

dLeft=diff(Left);
dLeft=(double(dLeft>0)-double(dLeft<0));
dRight=diff(Right);
dRight=(double(dRight>0)-double(dRight<0));
dTop=diff(Top);
dTop=(double(dTop>0)-double(dTop<0));
dBottom=diff(Bottom);
dBottom=(double(dBottom>0)-double(dBottom<0));

sdLeft=sum(abs(dLeft));
sdRight=sum(abs(dRight));
sdTop=sum(abs(dTop));
sdBottom=sum(abs(dBottom));

if(sdLeft>0)
    pixclr=Left(1);
    if(~pixclr)
        shift=1;
    end
    start=[find(Left~=Left(1),1,'first')-Left(1),1];
    if(sdTop>0)
        fin=[1,find(Top~=Top(1),1,'first')-Top(1)];
    elseif(sdRight>0)
        fin=[find(Right~=Right(1),1,'first')-Right(1),cols];
    elseif(sdBottom>0)
        fin=[rows,find(Bottom~=Bottom(1),1,'first')-Bottom(1)];
    else
        fin=[find(Left~=Left(1),1,'last')-Left(1),1];
    end
elseif(sdTop>0)
    start=[1,find(Top~=Top(1),1,'first')-Top(1)];
    if(sdRight>0)
        fin=[find(Right~=Right(1),1,'first')-Right(1),cols];
    elseif(sdBottom>0)
        fin=[rows,find(Bottom~=Bottom(1),1,'first')-Bottom(1)];
    elseif(sdLeft>0)
        fin=[find(Left~=Left(1),1,'first')-Left(1),1];
    else
        fin=[1,find(Top~=Top(1),1,'last')-Top(1)];
    end
elseif(sdRight>0)
    start=[find(Right~=Right(1),1,'first')-Right(1),cols];
    if(sdBottom>0)
        fin=[rows,find(Bottom~=Bottom(1),1,'first')-Bottom(1)];
    elseif(sdLeft>0)
        fin=[find(Left~=Left(1),1,'first')-Left(1),1];
    elseif(sdTop>0)
        fin=[1,find(Top~=Top(1),1,'first')-Top(1)];
    else
        fin=[find(Right~=Right(1),1,'last')-Right(1),cols];
    end        
elseif(sdBottom>0)
    start=[rows,find(Bottom~=Bottom(1),1,'first')-Bottom(1)];
    if(sdLeft>0)
        fin=[find(Left~=Left(1),1,'first')-Left(1),1];
    elseif(sdTop>0)
        fin=[1,find(Top~=Top(1),1,'first')-Top(1)];
    elseif(sdRight>0)
        fin=[cols,find(Right~=Right(1),1,'first')-Right(1)];
    else
        fin=[find(Bottom~=Bottom(1),1,'last')-Bottom(1),rows];
    end       
end

% Color=BW(1,1);% Color of first corner (first row, first col)
% for i=1:rows %starts loop row 1 and continues by increments of +one 
%     if (BW(i,1)~= Color) %if BW(i,1) does not equal Color   
%         point=[i-1,1];%then coordinates = point
%         break %loop breaks
%     end
% end
% if(i~=rows) % if color change did not happen in (1,1)
%     if isempty(start)%if start is empty
%         start= point;%start equals point(i,1)
%     else
%         fin=point; %otherwise fin equals point (i,1)
%     end
% end
% for i=2:cols
%     if (BW(1,i)~=Color)
%         point=[1,i-1];
%         break
%     end
% end
% if not(i==cols)
%     if isempty(start)
%         start= point;
%     else
%         fin=point;
%     end
% end
% Color=BW(rows,cols);
% for i=rows:-1:1
%     if (BW(i,cols)~= Color) 
%         point=[i-1,cols];
%         break
%     end
% end
% if not(i==1)
%     if isempty(start)
%         start= point;
%     else
%         fin=point;
%     end
% end
% for i=cols:-1:1
%     if (BW(rows,i)~=Color)
%         point=[rows,i-1];
%         break
%     end
% end
% if not(i==1)
%     if  isempty(start)
%         start= point;
%     else
%         fin=point;
%     end
% end
    %end endpts fxn
%     %------------------------------------------------------
    B=bwtraceboundary(BW,start,'S',8,Inf,'counterclockwise');%traces boundary from start in direction of steps
    if(norm(size(B))<1)
        fprintf(1,'Error with image: %s\n',pix.textdata{idx});
        continue;
    end
    X=abs(B(:,1)-fin(1));
    Y=abs(B(:,2)-fin(2));
    Z=X+Y;
    [loc,K]=find(Z==0);%index of fin point
    B=[B(1:loc,1),B(1:loc,2)];%only brd pts between start and fin
    num_r=2000;
    min_r = 2*sqrt(2);
    max_r = length(B(:,1))/(10*min_r);%2*length(B(:,1));%
    R=10.^linspace(log10(min_r),log10(max_r),2000);%defines step lengths
    
    matlab_loop=0;
    
    
    if(matlab_loop==0)
        %% run C code
        t1=tic;
        delete('boundary.dat','stepsizes.dat','params.dat','richdist.dat','polypts.dat');
        format long e;
        dlmwrite('boundary.dat',B,'precision', 16);
        dlmwrite('stepsizes.dat',R','precision', 16);
        fp=fopen('params.dat','w');
        fprintf(fp,'%d\n%d',length(B(:,1)),length(R));
        fclose(fp);
%         dlmwrite('params.dat',[length(B(:,1));length(R);]);
        base = pwd;
        t2=tic;
        system(['"',base,'/richDist"']);
        t2o=toc(t2);
        data_full=load('richdist.dat');        
        data=data_full(:,1:2);
        area=data_full(:,3);
        t1o=toc(t1);
        fprintf(1,'Full loop time: %.2e\n',t1o);
        fprintf(1,'C run time: %.2e\n',t2o);
    else
    %% Matlab loops
    Bi=B(:,1)+1i*B(:,2);%creates complex boundary pts
    k_len=1;
    steps_start=zeros(length(R),k_len);
    r_sizes=zeros(length(R),k_len);
    ave_steps=zeros(length(R),1);
    mode_steps=zeros(length(R),1);
    parfor (j=1:length(R)) %loop through all step sizes
    r=R(j);
    r_size=zeros(k_len,1);
    tot_steps=zeros(k_len,1);
    for k = 1:k_len; %loops through all starting pts
        step_max=5*r;
        tot_steps(k)=richardsonDistance(B(:,1),B(:,2),r,step_max); %#ok<PFBNS>
        r_size(k)=r;
    end
    steps_start(j,:)=tot_steps;
    r_sizes(j,:)=r_size;
    ave_steps(j)=mean(tot_steps);
    mode_steps(j)=mode(tot_steps);

    end
    data = [r_sizes,ave_steps];
%tot_realz=[tot_real./pc];
    end
tend=toc(tstart);
fprintf(1,'Runtime (total): %.2e\n',tend);

%% Post processing
% data=data(~isnan(data(:,2)),:);
log_real_stepsize = log10(data(:,1)/pc);
log_real_lengths = log10(data(:,2)/pc);
p=polyfit(log_real_stepsize,log_real_lengths,1);
display(p);
psteps=polyval(p,log_real_stepsize);
err=psteps-log_real_lengths;
SSE=err'*err;
Rsquared=1-SSE/length(err);

% plotting
plot(log_real_stepsize,log_real_lengths,'.');
hold on;
set(0,'defaulttextinterpreter','latex');
plot(log_real_stepsize,psteps,'m','LineWidth',3);
set(gca,'FontSize',15);
ylabel('Perimeter: $\log_{10}$(cm)');
xlabel('Step Size: $\log_{10}$(cm)');
drawnow;
fh=fopen([pix.textdata{idx},'.txt'],'w');
fprintf(fh,'1-D=%.15e, log(M)=%.15e\n',p(1),p(2));
fclose(fh);
fh=fopen([pix.textdata{idx},'.csv'],'w');
fprintf(fh,'Ruler in cm, Perimeter in cm, Area in cm2\n');
for i=1:length(data(:,1))
    fprintf(fh,'%.15e, %.15e, %.15e\n',data(i,1)/pc,data(i,2)/pc,area(i)/pc/pc);
end
fclose(fh);


% legend({'data';sprintf('y=%.3fx+%.3f\n$R^2=$%.3f',...
%     p(1),p(2),1-SSE/length(err))},'Interpreter','Latex')

% pts = load('polypts.dat');
% figure();
% plot(B(:,1),B(:,2),'b.-')
% hold on
% plot(pts(:,1),pts(:,2),'r.-')
% axis equal
end
