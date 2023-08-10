close all;
clear;
clc;
dfile='random_deployment.txt';
if exist(dfile,'file')
    delete(dfile);
end
diary (dfile)
diary on
%%%%%%%%%%%%%%%%%%%% Network Establishment Parameters %%%%%%%%%%%%%%%%%%%%
% Field Dimensions in meters %
%Reading network parameters from the config text file
fid = fopen('config.txt');
C = textscan(fid, '%[^= ]%*[= ]%f', 'CommentStyle', '%');
fclose(fid);
xm = C{2}(strcmp(C{1}, 'xm')); %length of the network field in x-axis
ym = C{2}(strcmp(C{1}, 'ym')); %length of the network field in y-axis
r = C{2}(strcmp(C{1}, 'r')); %transmission radius of sensor nodes
nd = C{2}(strcmp(C{1}, 'nd')); % node density of network region
packets = C{2}(strcmp(C{1}, 'packets')); %number of packets to be transmitted by the source node
%Parameter settings%
x=0;
y=0; % added for better display results of the plot
% threshold distance  %
d0=87;
%number of nodes%
n=xm*xm*nd;
% Coordinates of the Sink (location is predetermined in this simulation) %
sinkx=xm/2;
sinky=ym/2;

%%% Energy Values %%%
% Initial Energy of a Node (in Joules) % 
Eo=0.5; % units in Joules
% Energy required to run circuity (both for transmitter and receiver) %
Eelec=50*10^(-9); % units in Joules/bit
Efs=10*10^(-12);
% Transmit Amplifier Types %
Eamp=0.0013*10^(-12); % units in Joules/bit/m^2 (amount of energy spent by the amplifier to transmit the bits)
%average energy consumption
aec=0;
% Size of data package %
b=1024; % units in bits
% Suggested percentage of cluster head %
p=0.05; % a 5 percent of the total amount of nodes used in the network is proposed to give good results 
%%%%%%%%%%%%%%%%%%%%%%%%%%% End of Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf ( '\n Value for x-axis is %d', xm);
fprintf ( '\n Value for y-axis is:  %d', ym);
fprintf ( '\n The grid size is:  %d', r);
fprintf ( '\n The network density is:  %d', nd);
            %%% Creation of the Wireless Sensor Network %%%
% Plotting the WSN %
for i=1:n
    SN(i).id=i;	% sensor's ID number
    SN(i).x=rand(1,1)*xm;	% X-axis coordinates of sensor node
    SN(i).y=rand(1,1)*ym;	% Y-axis coordinates of sensor node
    SN(i).E=Eo;     % nodes energy levels (initially set to be equal to "Eo"
    SN(i).role=0;   % node acts as normal if the value is '0', if elected as a cluster head it  gets the value '1' (initially all nodes are normal)
    SN(i).cluster=0;	% the cluster which a node belongs to
    SN(i).cond=1;	% States the current condition of the node. when the node is operational its value is =1 and when dead =0
    SN(i).rop=0;	% number of rounds node was operational
    SN(i).rleft=0;  % rounds left for node to become available for Cluster Head election
    SN(i).dtch=0;	% nodes distance from the cluster head of the cluster in which he belongs
    SN(i).dts=0;    % nodes distance from the sink
    SN(i).tel=0;	% states how   many times the node was elected as a Cluster Head
    SN(i).rn=0;     % round node got elected as cluster head
    SN(i).chid=0;   % node ID of the cluster head which the "i" normal node belongs to
    SN(i).packets=60; %maximum number of packets that can be sent in the transmission
    SN(i).cttrsmn=0; %counter for transmission is for counting how many times the node is selected for transmission
    hold on;
    figure(1)
    plot(x,y,xm,ym,SN(i).x,SN(i).y,'ob',sinkx,sinky,'*r');
    text(sinkx, sinky,'sink', 'Color','b');
    title 'Wireless Sensor Network';
    xlabel '(m)';
    ylabel '(m)';
    grid on
    xd=40;
  tickvalues= min(x):xd:max(xm);
    set(gca, 'XTick',tickvalues);
    yd=40;
    tickValues2 = min(y):yd:max(ym);
    set(gca,'YTick',tickValues2); 
end
saveas(gcf,'randomdeployment.jpg');
xvector = 0:r:xm;
yvector = 0:r:ym;
k=r;
for i=1:n
    
    fprintf('\n for node=%d', i);
    point= [ SN(i).x,SN(i).y];
fprintf( '\n coordinates are= %d, %d ', SN(i).x, SN(i).y);
x_grid = find( xvector <= point(1), 1, 'last');
y_grid = find(yvector<= point(2), 1, 'last');
x_grid=(x_grid*k)-k;
y_grid=(y_grid*k);

SN(i).x_grid = [x_grid x_grid+k];
SN(i).y_grid = [y_grid y_grid-k];
 SN(i).upgrid= [x_grid  y_grid ];
 SN(i).lpgrid = [x_grid+k y_grid-k];
fprintf('Grid value:(%d, %d),(%d, %d) ', SN(i).x_grid(1), SN(i).y_grid(1),SN(i).x_grid(2), SN(i).y_grid(2));
end
nofrings= sinkx/k;
fprintf('\n No. of rings=%d', nofrings);
totalcells=0;total=0;
for i=1:nofrings
    R(i).no=i;
    upx=sinkx-(k*i);
    upy= sinky+(k*i);
    R(i).up=[upx, upy];
    lpx=sinkx+(k*i);
    lpy=sinky-(k*i);
    R(i).lp=[lpx ,lpy];
   fprintf('\n \n');
    fprintf('\n For %dth ring',R(i).no);
    fprintf('\n Upper coordinates are= %d, %d', R(i).up(1),R(i).up(2));
    fprintf('\n Lower coordinates are= %d ,%d', R(i).lp(1),R(i).lp(2));
    R(i).nocx= (lpx-upx)/k;            % no.of cells along x axis
    R(i).nocy= (upy-lpy)/k;            %no. of cells along y axis
    R(i).tcells= R(i).nocx*R(i).nocy;     %to calculate total cells within that specified upper and lower coords
   if i==1
    R(i).totalcells= R(i).tcells;  
   else 
       R(i).totalcells= R(i).tcells-R(i-1).tcells;  % to calculate total cells within that ring leaving out the cells in inner ring
   end
       
    fprintf('\n no. of cells= %d', R(i).totalcells);
t=1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%  finding cells within each ring along with their coords
    %%%%%%%%%%   C is a struct containing the cords of each cell in their
    %%%%%%%%%%   specific ring.
    fprintf(' \n The cells along with their upper and lower coordinates are: \n'); 
    
    for t=1: R(i).totalcells   %for  each cell inside the ring
        if(t==1)
            cux1=R(i).up(1);
            cuy1=R(i).up(2);
            R(i).C(t).up= [cux1 cuy1];
            clx2= R(i).up(1)+k;
            cly2= R(i).up(2)-k;
            R(i).C(t).lp=[clx2 cly2];
            R(i).C(t).cno=1;
            
        elseif (clx2< R(i).lp(1)) && (cuy1==R(i).up(2))
                 cux1=R(i).C(t-1).up(1)+k;
                cuy1=R(i).C(t-1).up(2);
                R(i).C(t).up= [cux1 cuy1];
                clx2= R(i).C(t-1).lp(1)+k;
                cly2= R(i).C(t-1).lp(2);
                 R(i).C(t).lp=[clx2 cly2];
                 R(i).C(t).cno=R(i).C(t-1).cno+1;
        elseif (cly2> R(i).lp(2)) && (clx2==R(i).lp(1))
                 cux1=R(i).C(t-1).up(1);
                cuy1=R(i).C(t-1).up(2)-k;
                R(i).C(t).up= [cux1 cuy1];
                clx2= R(i).C(t-1).lp(1);
                cly2= R(i).C(t-1).lp(2)-k;
                 R(i).C(t).lp=[clx2 cly2];
                 R(i).C(t).cno=R(i).C(t-1).cno+1;
        elseif (cux1 > R(i).up(1)) && (cly2==R(i).lp(2))
               cux1=R(i).C(t-1).up(1)-k;
                cuy1=R(i).C(t-1).up(2);
                R(i).C(t).up= [cux1 cuy1];
                clx2= R(i).C(t-1).lp(1)-k;
                cly2= R(i).C(t-1).lp(2);
                 R(i).C(t).lp=[clx2 cly2];
                 R(i).C(t).cno=R(i).C(t-1).cno+1;
         else
                 cux1=R(i).C(t-1).up(1);
                cuy1=R(i).C(t-1).up(2)+k;
                R(i).C(t).up= [cux1 cuy1];
                clx2= R(i).C(t-1).lp(1);
                cly2= R(i).C(t-1).lp(2)+k;
                 R(i).C(t).lp=[clx2 cly2];
                 R(i).C(t).cno=R(i).C(t-1).cno+1;
        end
        fprintf('  %dth:(%d, %d)  & (%d,%d),',R(i).C(t).cno, R(i).C(t).up(1),R(i).C(t).up(2),R(i).C(t).lp(1), R(i).C(t).lp(2));
          %%%%%%%%%to calculate the nodes within each cell%%%%%%%%%%%%%%% 
        nod=0;           %to count number of nodes inside each cell
        R(i).C(t).nodeloc=[];R(i).C(t).nodeid=[];
      for  s=1:n 
          if (SN(s).x>=R(i).C(t).up(1)) && (SN(s).x<=R(i).C(t).lp(1)) && (SN(s).y<=R(i).C(t).up(2)) && (SN(s).y>=R(i).C(t).lp(2))
          nod=nod+1; 
           R(i).C(t).nodeloc=[ R(i).C(t).nodeloc SN(s).x  SN(s).y]; % R is a struct and nodeloc stores the location of each node  within the ring 
            R(i).C(t).nodeid= [ R(i).C(t).nodeid SN(s).id];
       
          end
          R(i).C(t).nodes=nod;
      end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%calculating number of nodes within each ring and their location
    
    R(i).nodeloc=[];R(i).nodeid=[];
    count=0; %to calculate the no. of nodes within a ring
        for j=1:n                    %loop to find out nodes inside the ring
        if(i==1)
            if(SN(j).x>=R(i).up(1) & SN(j).x<=R(i).lp(1) & SN(j).y <=R(i).up(2) & SN(j).y>=R(i).lp(2))
            count=count+1;
             SN(j).rno=i;
            R(i).nodeloc=[R(i).nodeloc SN(j).x  SN(j).y]; % R is a struct and nodeloc stores the location of each node  within the ring 
            R(i).nodeid= [R(i).nodeid SN(j).id];  
%             fprintf('%dth:(%d, %d) , ', SN(j).id, R(count).nodeloc(1), R(count).nodeloc(2));
            end
        else
           if((SN(j).x>R(i).up(1) && SN(j).x<R(i).lp(1) & SN(j).y <R(i).up(2) & SN(j).y>R(i).lp(2)) & ~(SN(j).x>R(i-1).up(1) & SN(j).x<R(i-1).lp(1) & SN(j).y <R(i-1).up(2) & SN(j).y>R(i-1).lp(2)))
             count=count+1;
              SN(j).rno=i;
              R(i).nodeloc=[R(i).nodeloc SN(j).x  SN(j).y];
                R(i).nodeid= [R(i).nodeid SN(j).id];
%               fprintf('%dth :(%d, %d) ', SN(j).id, R(count).nodeloc(1), R(count).nodeloc(2));
           end
          
        end     
        end
    R(i).ct=count;         %calculate no. of nodes within each ring
    fprintf('\n total no. of nodes in this ring is %d  ', R(i).ct);
    total=R(i).ct+total;
    fprintf('\n the nodes in this ring with id and loc are :\n' );
    for cn=1:count
            fprintf('%dth:(%d, %d) , ', R(i).nodeid(cn), R(i).nodeloc(cn), R(i).nodeloc(cn+1));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%finding neighbor nodes of each cell in each ring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     maxn=8; 
      for j=1: R(i).totalcells
       ncno=0;
     fprintf('\n for cell %d the neighbor cells are: ', R(i).C(j).cno);
      for nn=1: maxn
              if (nn==1)
                  ncupx= R(i).C(j).up(1)-k;
                  ncupy=R(i).C(j).up(2)+k;
                 R(i).C(j).NC(nn).up=[ncupx ncupy];
                  nclpx=R(i).C(j).up(1);
                  nclpy=R(i).C(j).up(2);
                   R(i).C(j).NC(nn).lp=[nclpx  nclpy];
                    R(i).C(j).NC(nn).ncno=ncno+1;
                  
              elseif (ncupx<(R(i).C(j).up(1)+k)) && (nclpy==R(i).C(j).up(2))
                  ncupx= R(i).C(j).NC(nn-1).up(1)+k;
                  ncupy=  R(i).C(j).NC(nn-1).up(2);
                   R(i).C(j).NC(nn).up= [ncupx ncupy];
                  nclpx=  R(i).C(j).NC(nn-1).lp(1)+k;
                  nclpy= R(i).C(j).NC(nn-1).lp(2);
                   R(i).C(j).NC(nn).lp=[nclpx nclpy];
                   R(i).C(j).NC(nn).ncno= R(i).C(j).NC(nn-1).ncno+1;
                  
              elseif(nclpy > (R(i).C(j).lp(2)-k ) && (nclpx== (R(i).C(j).lp(1)+k)))
                  ncupx= R(i).C(j).NC(nn-1).up(1);
                  ncupy=  R(i).C(j).NC(nn-1).up(2)-k;
                   R(i).C(j).NC(nn).up= [ncupx ncupy];
                  nclpx= R(i). C(j).NC(nn-1).lp(1);
                  nclpy= R(i).C(j).NC(nn-1).lp(2)-k;
                   R(i).C(j).NC(nn).lp=[nclpx nclpy];
                   R(i).C(j).NC(nn).ncno=R(i).C(j).NC(nn-1).ncno+1;
                  
              elseif (ncupx>(R(i).C(j).up(1)-k))
                   ncupx=R(i). C(j).NC(nn-1).up(1)-k;
                  ncupy=  R(i).C(j).NC(nn-1).up(2);
                   R(i).C(j).NC(nn).up= [ncupx ncupy];
                  nclpx= R(i). C(j).NC(nn-1).lp(1)-k;
                  nclpy= R(i).C(j).NC(nn-1).lp(2);
                   R(i).C(j).NC(nn).lp=[nclpx nclpy];
                   R(i).C(j).NC(nn).ncno=R(i).C(j).NC(nn-1).ncno+1;
                  
              else
                  ncupx= R(i).C(j).NC(nn-1).up(1);
                  ncupy=  R(i).C(j).NC(nn-1).up(2)+k;
                   R(i).C(j).NC(nn).up= [ncupx ncupy];
                  nclpx=  R(i).C(j).NC(nn-1).lp(1);
                  nclpy= R(i).C(j).NC(nn-1).lp(2)+k;
                   R(i).C(j).NC(nn).lp=[nclpx nclpy];
                   R(i).C(j).NC(nn).ncno=R(i).C(j).NC(nn-1).ncno+1;
              end
              fprintf('(%d, %d)(%d, %d)&',  R(i).C(j). NC(nn).up(1),  R(i).C(j).NC(nn).up(2) ,  R(i).C(j). NC(nn).lp(1),   R(i).C(j).NC(nn).lp(2)); 
%               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       %%%%%%%%%%%%%%%%%%%to find nodes  present in the neighbor cells%%%%%%%
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       for ss=1:n
%           if (R(i).C(j).NC(nn).up(1)==SN(ss).upgrid(1)) & (R(i).C(j).NC(nn).up(2)==SN(ss).upgrid(2))
%               R(i).C(j).NC(nn).nodeid =SN(ss).id; 
%           else
%               R(i).C(j).NC(nn).nodeid =0;
%           end
%       end
      end
      end
end
 %%%%for each ring find the clusterhead in each cell
 clheads=0;
    for r=1:nofrings  %for each ring
        for c=1:R(r).totalcells  %for each cell in ring
            if (R(r).C(c).nodes >1)
                R(r).C(c).chead=R(r).C(c).nodeid(1);
                for no=2:R(r).C(c).nodes   %for comparing nodes present in cell
                    if (SN(R(r).C(c).nodeid(no)).E) > (SN(R(r).C(c).chead).E)
                            R(r).C(c).chead=R(r).C(c).nodeid(no);  
                    else if (SN(R(r).C(c).nodeid(no)).E) < (SN(R(r).C(c).chead).E)
                             R(r).C(c).chead=R(r).C(c).chead;    
                        else 
                            R(r).C(c).chead=R(r).C(c).nodeid(no);      
                        end
                    end
                end
                clheads=clheads+1;
                SN(R(r).C(c).chead).role=1;
            elseif (R(r).C(c).nodes ==1)  %if there is only one node, that will be the clhead
              R(r).C(c).chead=R(r).C(c).nodeid(1); 
              clheads=clheads+1;
              SN(R(r).C(c).nodeid(1)).role=1;
            else 
                R(r).C(c).chead=0; 
            end
        end
        
    end
    
   %%%%%%%%%%%%%%%%%to find the clusterhead in each neighboring node%%%%
     for i=1:nofrings
      for j=1: R(i).totalcells
      for nn=1:maxn
      if (i~=nofrings)
          for u=1:R(i+1).totalcells
              if ((R(i).C(j). NC(nn).up(1)== R(i+1).C(u).up(1)) & (R(i).C(j). NC(nn).up(2)== R(i+1).C(u).up(2)))
                 R(i).C(j).NC(nn).nodeid= R(i+1).C(u).nodeid;
                 R(i).C(j).NC(nn).nodes= R(i+1).C(u).nodes;
                  R(i).C(j).NC(nn).chead= R(i+1).C(u).chead;
                   R(i).C(j).NC(nn).rno=i+1;
                    R(i).C(j).NC(nn).cno= R(i+1).C(u).cno;
              end
          end
      end
      for v=1:R(i).totalcells
          if ((R(i).C(j). NC(nn).up(1)== R(i).C(v).up(1)) & (R(i).C(j). NC(nn).up(2)== R(i).C(v).up(2)))
              R(i).C(j).NC(nn).nodeid= R(i).C(v).nodeid;
              R(i).C(j).NC(nn).nodes= R(i).C(v).nodes;
              R(i).C(j).NC(nn).chead= R(i).C(v).chead;
                  R(i).C(j).NC(nn).rno=i;
                  R(i).C(j).NC(nn).cno= R(i).C(v).cno;
          end
      end
      if (i~=1)
      for z=1:R(i-1).totalcells
          if ((R(i).C(j). NC(nn).up(1)== R(i-1).C(z).up(1)) & (R(i).C(j). NC(nn).up(2)== R(i-1).C(z).up(2)))
              R(i).C(j).NC(nn).nodeid= R(i-1).C(z).nodeid;
              R(i).C(j).NC(nn).nodes= R(i-1).C(z).nodes;
              R(i).C(j).NC(nn).chead= R(i-1).C(z).chead;
                  R(i).C(j).NC(nn).rno=i-1;
                  R(i).C(j).NC(nn).cno= R(i-1).C(z).cno;
          end
      end
      end
      end 
    if i==nofrings
          for nk=1:maxn
           t= isempty(R(i).C(j).NC(nk).chead);
           if t==1
                R(i).C(j).NC(nk).nodes=0;
               R(i).C(j).NC(nk).nodeid=0;
               R(i).C(j).NC(nk).chead=0;
               R(i).C(j).NC(nk).rno=0;
               R(i).C(j).NC(nk).cno=0;
           end
          end
    end
      end
     end

fprintf('\n total no. of nodes %d \n', total);
diary off
