close all;
clear;
clc;
dfile='shortest_path_routing.txt';
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
% Number of Nodes in the field %
init_n=0; %to start the numbering of nodes
% n=480;
% Number of Dead Nodes in the beggining %
dead_nodes=0;
% threshold distance  %
d0=87;
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
xd=r;%Grid size in x-axis
yd=r; %Grid size in y-axis
k=r;
%%%%%%%%%%%%%%%%%%%%%%%%%%% End of Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf ( '\n Value for x-axis is %d', xm);
fprintf ( '\n Value for y-axis is:  %d', ym);
fprintf ( '\n The grid size is:  %d', r);
fprintf ( '\n The network density is:  %d', nd);
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
    end
    % to decide about the color of sensor nodes in the rings %%%%%%%%%%
    if mod(i,2)==0
        col='ob';
    else
        col='og';
    end
        % Plotting the WSN %
    nc=round(xd*yd*nd); %find out number of nodes in each cell
    for tn=1: R(i).totalcells %to plot sensor nodes in each cell
        for z=init_n+1:init_n+nc   % for each sensor node that gets added by every 5th+1 sensor in each cell
            SN(z).id=z;	% sensor's ID number
            SN(z).x=randi([R(i).C(tn).up(1)+1,R(i).C(tn).lp(1)-1],1,1); % rand(1,1)*xd;	% X-axis coordinates of sensor node, + and - to not allow any node on the boundary
            SN(z).y=randi([R(i).C(tn).lp(2)+1,R(i).C(tn).up(2)-1],1,1);%rand(1,1)*yd;	% Y-axis coordinates of sensor node
            SN(z).E=Eo;     % nodes energy levels (initially set to be equal to "Eo"
            SN(z).role=0;   % node acts as normal if the value is '0', if elected as a cluster head it  gets the value '1' (initially all nodes are normal)
            SN(z).cluster=0;	% the cluster which a node belongs to
            SN(z).cond=1;	% States the current condition of the node. when the node is operational its value is =1 and when dead =0
            SN(z).rop=0;	% number of rounds node was operational
            SN(z).rleft=0;  % rounds left for node to become available for Cluster Head election
            SN(z).dtch=0;	% nodes distance from the cluster head of the cluster in which he belongs
            SN(z).dts=0;    % nodes distance from the sink
            SN(z).tel=0;	% states how   many times the node was elected as a Cluster Head
            SN(z).rn=0;     % round node got elected as cluster head
            SN(z).chid=0;   % node ID of the cluster head which the "i" normal node belongs to
            SN(z).ct_pkt=0;% for counting the no.of times a node forwards a packet
            SN(z).cttrsmn=0; %counter for transmission is for counting how many times the node is selected for transmission
            init_n=z; %to start the count of nodes in thr nrxt cell from next sensor number and not from 1 again
%             hold on;
%             figure(1)
%             plot(x,y,xm,ym,SN(z).x,SN(z).y,col,sinkx,sinky,'*r');
%             text(sinkx, sinky,'sink', 'Color','b');
%             title 'Wireless Sensor Network'
%             xlabel '(m)';
%             ylabel '(m)';
%              grid on
%              
        end
        %     xd=40;
%         tickvalues= min(x):xd:max(xm);
%         set(gca, 'XTick',tickvalues,'GridLineStyle','--');
% %         set(gca,'GridLineStyle','--')
%         %     yd=40
%         tickValues2 = min(y):yd:max(ym);
%         set(gca,'YTick',tickValues2,'GridLineStyle','--');
%         
%         xtickangle(45);
% %         ytickangle(45);
    end
end
n=z; %no.of nodes  
        
        %%%%%%%%%to calculate the nodes within each cell%%%%%%%%%%%%%%% 
        for i=1:nofrings
            for t=1:R(i).totalcells
                 R(i).C(t).nodeloc=[];
                 R(i).C(t).nodeid=[];
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
         end
    
    %%%calculating number of nodes within each ring and their location
    for i=1:nofrings
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
    fprintf('\n \n total no. of nodes in %d ring is %d  ', i, R(i).ct);
    total=R(i).ct+total;
    fprintf('\n the nodes in this ring with id and loc are :\n' );
    for cn=1:count
            fprintf('%dth:(%d, %d) , ', R(i).nodeid(cn), R(i).nodeloc(cn), R(i).nodeloc(cn+1));
    end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%finding neighbor cells of each cell in each ring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nofrings
%    fprintf('\n for ring %d : ', i);     
maxn=8; 
      for j=1: R(i).totalcells
       ncno=0;
%      fprintf('\n for cell %d the neighbor cells/grids are: ', R(i).C(j).cno);
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
%               fprintf('(%d, %d)(%d, %d)&',  R(i).C(j). NC(nn).up(1),  R(i).C(j).NC(nn).up(2) ,  R(i).C(j). NC(nn).lp(1),   R(i).C(j).NC(nn).lp(2)); 
      end
      end
end
fprintf('\n total no. of nodes %d \n', total);


%%%%%%%%%%%%%%%%%  Hop count for each node   %%%%%%%%%%%%%%%%%55
for i=1:n
    
        SN(i).E=Eo;
        SN(i).hc= SN(i).rno;
end

     %%%%%%%%%%%%%%%randomly generating source node%%%%%%%%%%%%%%%%
     src=randi([1,n],1,1);
     fprintf(" \n Source node= %d", src);
     fprintf('  Source node location:(%d, %d) ', SN(src).x, SN(src).y);
     fprintf(" \n Ring no.= %d", SN(src).rno);
     for cell=1: R(SN(src).rno).totalcells
         for nod=1:R(SN(src).rno).C(cell).nodes
             if (R(SN(src).rno).C(cell).nodeid(nod)==src)
                 sncn=R(SN(src).rno).C(cell).cno; % sncn is source node cell no. to store the cell no. of source node
                 fprintf(" \n cell no.= %d",sncn);
             end
         end
     end 
     rnd=packets; % total number of rounds
     ct_assetcaptured=0;% counter to count the number of times asset captured
     assetcapturedround=0;
    for round=1:rnd
          %   round
         fprintf("\n Sending %d th packet ",round);
         clheads=0;
         total_hops=0; %to find the total number of hops per packet 
         tec=0; %to find the total energy consumption per packet
        
    %%%%for each ring find the clusterhead in each cell
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
          if ((R(i).C(j). NC(nn).up(1)== R(i-1).C(z).up(1)) && (R(i).C(j). NC(nn).up(2)== R(i-1).C(z).up(2)))
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
  
%%%%Shortest path routing%%%%%%
pnrn= SN(src).rno;
pncn=sncn;
pn=R(pnrn).C(pncn).chead;
Rd(round).SPR(1).nodeid=pn;  %%%to store the source node as the first starting node in the structure
Rd(round).SPR(1).cell=pncn;
Rd(round).SPR(1).ring=SN(pn).rno; 
ct_spr=1;
while(pnrn>1)
    nbrs=0;
    for nn=1:maxn
        if (R(pnrn).C(pncn).NC(nn).rno==pnrn-1)
            nbrs=nbrs+1;
            Neighborlist(nbrs).id=ct_spr;
            Neighborlist(nbrs).pn=R(pnrn).C(pncn).NC(nn).chead;
            Neighborlist(nbrs).rno=R(pnrn).C(pncn).NC(nn).rno;
           Neighborlist(nbrs).cell=R(pnrn).C(pncn).NC(nn).cno;
        end
    end       
pn_random=randi([1,nbrs],1,1);
next_pn=Neighborlist(pn_random).pn; %next_pn denotes the next primary node that will receive the packet
cordinates=[SN(pn).x, SN(pn).y; SN(next_pn).x,SN(next_pn).y]; %stores the two cordinates we have to compute the distance
d=pdist(cordinates,'euclidean');
if d<d0 
    et= Eelec*b+ Efs*b*(d^2);
else
    et=Eelec*b+Efs*b*(d^4);
end
er=Eelec*b;
SN(next_pn).E=SN(next_pn).E-(et+er); 
tec=tec+et+er;
pn=next_pn;
pncn=Neighborlist(pn_random).cell;
pnrn=Neighborlist(pn_random).rno;
SN(pn).ct_pkt=SN(pn).ct_pkt+1;
ct_spr=ct_spr+1;
Rd(round).SPR(ct_spr).nodeid=pn;
Rd(round).SPR(ct_spr).cell=pncn;
Rd(round).SPR(ct_spr).ring=SN(pn).rno; 
end
%%%%%%%%%%%recording hop count  for shortest path routing %%%%%%%%%%%%%%%%%%%%%%
Rd(round).hop_spr=ct_spr;
%%%%%%%%%%recording total hops traversed by each packet%%%%%%%%%%%%%
total_hops=ct_spr;
Rd(round).total_hops=total_hops;
%%%%%%%%%recording total energy consumed in the entire process of routing%
Rd(round).total_energy=tec;
    
%%%%%%%%%%%%%%%%%%  ADVERSARIAL CODE %%%%%%%%%%%%%%%%%%%%
for pc=1:total_hops%%to assign all the nodes in a single structure, hops+1 as first node is the src node and counter for spr starts from 1
    Rd(round).Allnodes(pc).pn_nid=Rd(round).SPR(pc).nodeid;
    Rd(round).Allnodes(pc).ring=Rd(round).SPR(pc).ring; 
    Rd(round).Allnodes(pc).cell=Rd(round).SPR(pc).cell;
end
       
ctr_advbcktrk=0; %%%%counter for choice of no. of adversary backtrack moves
if round==1
    adv=Rd(round).Allnodes(pc).pn_nid; %%for 1st packet, adv will start backtracking from the last node that traverses the packet
    adrn=SN(adv).rno;
    adcn=Rd(round).Allnodes(pc).cell;
    ctbtk=1; %counter for counting no. of backtracks
    adversarylocation(ctbtk).node=adv;
    adversarylocation(ctbtk).ring=adrn;
    adversarylocation(ctbtk).cell=adcn;
    adversarylocation(ctbtk).round=round;
    %btk=ctp3-1; %btk is for backtracking nodes counter
else
    for s=total_hops:-1:2
              if (Rd(round).Allnodes(s).ring==adrn) & (Rd(round).Allnodes(s).cell==adcn)
                  adv=Rd(round).Allnodes(s-1).pn_nid;
                  adrn=Rd(round).Allnodes(s-1).ring;
                  adcn=Rd(round).Allnodes(s-1).cell;
                  ctbtk=ctbtk+1;
                  adversarylocation(ctbtk).node=adv;
                  adversarylocation(ctbtk).ring=adrn;
                  adversarylocation(ctbtk).cell=adcn;
                  adversarylocation(ctbtk).round=round;
                  break;
              end
    end
end
if (adrn==SN(src).rno) && (adcn==sncn)
    ct_assetcaptured=ct_assetcaptured+1;
    assetcapturedround(ct_assetcaptured)=round; 
    adv=SN(src).id;
    adrn=SN(src).rno;
    adcn=sncn;
    adversarylocation(ctbtk).node=adv;
    adversarylocation(ctbtk).ring=adrn;
    adversarylocation(ctbtk).cell=adcn;
    adversarybacktrack=ctbtk; %%%to count no. of backtracks done by adversary
    safetyperiod=assetcapturedround;
    break;
end


    end
    total_rounds=round;
    %%calculate entropy after all the rounds for each simulation
     entropy=0;
     for sensor=1:n
         SN(sensor).pi=SN(sensor).ct_pkt/ total_rounds;
         SN(sensor).ratio=SN(sensor).pi*log2(1/SN(sensor).pi);
         SN(sensor).ratio(isnan(SN(sensor).ratio))=0;%%to convert NaN values to 0otherwise the addition wouldnt happen
         entropy=entropy+(SN(sensor).ratio);
     end
     %%%computing average transmission delay and total energy consumption%
     total_latency=0;
     transmission_latency=0;
     total_energy=0;
 for r=1:total_rounds
     total_latency=Rd(r).total_hops+total_latency;
     total_energy=Rd(r).total_energy+total_energy;
 end
 transmission_latency=total_latency/total_rounds;
 fprintf('\nSafety period at the end of  the simulation is  : %d',safetyperiod);
 fprintf('\nEntropy at the end of the simulation is  : %d',entropy);
 fprintf('\nTransmission latency at the end of the simulation is  : %d',transmission_latency);
 fprintf('\nTotal Energy consumption at the end of  the simulation is  : %d', total_energy);

diary off 
 
 
     
     
         
    
     
