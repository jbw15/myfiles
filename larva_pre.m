
clear all
clc
dt = 5e-4; t = dt:dt:1;num = length(t);
ds = 0.01; s = 0:ds:1; nl_totl = length(s) - 1; ms = ds/2:ds:1-ds/2;
dm = ds*ones(1,nl_totl); m = sum(dm);
omega = 2*pi; lambda =5; k = 2*pi/lambda;
% a1 = -1; a2 = 1; a3 = -0.2;
% am = a1*s.^2 + a2*s +a3;
a = 1; am = a*(s-2).^2+3;
nw = 6;
for nt = 1:num
    t = nt*dt;
%     if t<0.5
       curv(nt,:) = am.*sin(k*s-omega*t);
%     else
%        curv(nt,:) = am.*sin(k*s-omega*t); 
%     end
%     figure
%     plot(s,curv(nt,:))
    for nl = 2:nl_totl+1
         xx = [s(1:nl)];the = [curv(nt,1:nl)];
        theta(nt,nl) = trapz(xx,the)+pi;
         sxx = cos(theta(nt,1:nl));syy = sin(theta(nt,1:nl));
        x0(nt,nl)=trapz(xx,sxx);
        y0(nt,nl)=trapz(xx,syy);
    end
%     figure
%     plot(x0(nt,:),y0(nt,:))
%     hold on
%     plot1(x0(nt,1),y0(nt,1),0.005)
    X0(nt,:)=0.5*(x0(nt,1:nl_totl)+x0(nt,2:nl_totl+1));
    Y0(nt,:)=0.5*(y0(nt,1:nl_totl)+y0(nt,2:nl_totl+1));    
    cx(nt) = sum(dm.*X0(nt,:))/m;
    cy(nt) = sum(dm.*Y0(nt,:))/m;
    
    x1(nt,:)=x0(nt,:)-cx(nt);
    y1(nt,:)=y0(nt,:)-cy(nt);
    X1(nt,:)=X0(nt,:)-cx(nt);
    Y1(nt,:)=Y0(nt,:)-cy(nt);
%     figure
%     plot(x1(nt,:),y1(nt,:))
%     hold on
%     plot1(x1(nt,1),y1(nt,1),0.005)
%     axis([-0.5 0.5 -0.5 0.5])
%     i=num2str(nt);
%     saveas(gcf,[i,'.png']);
    cx1(nt)=sum(dm.*X1(nt,:))/m;
    cy1(nt)=sum(dm.*Y1(nt,:))/m;

end
for nt = 1:num
    X1(nt,:)=X1(nt,:)-cx1(nt);
    Y1(nt,:)=Y1(nt,:)-cy1(nt);
    x1(nt,:)=x1(nt,:)-cx1(nt);
    y1(nt,:)=y1(nt,:)-cy1(nt);
    ic1(nt,:)=dm.*(X1(nt,:).^2+Y1(nt,:).^2);
    Icom1(nt)=sum(ic1(nt,:));
    r_2(nt,:) = X1(nt,:).^2 + Y1(nt,:).^2;
    r(nt,:) = sqrt(r_2(nt,:));   
end
% for nt = 2:num-1
%     vx(nt,:) = (X1(nt+1,:) - X1(nt-1,:))/(2*dt);
%     vy(nt,:) = (Y1(nt+1,:) - Y1(nt-1,:))/(2*dt);
%     vr(nt,:) = (r(nt+1,:) - r(nt-1,:))/(2*dt);
% end
% nt = 1;
% vx(nt,:) = 2*vx(nt+1,:) - vx(nt+2,:); vy(nt,:) = 2*vy(nt+1,:)-vy(nt+2,:);vr(nt,:) = 2*vr(nt+1,:)-vr(nt+2,:);
% nt = num;
% vx(nt,:) = 2*vx(nt-1,:) - vx(nt-2,:); vy(nt,:) = 2*vy(nt-1,:)-vy(nt-2,:);vr(nt,:) = 2*vr(nt-1,:)-vr(nt-2,:);
% for nt = 1:num
%     for nl = 1:nl_totl
%         r1(nl,:) = [X1(nt,nl),Y1(nt,nl),0];
%         LA(nl) = norm(r1(nl,:));
%         if nt<num
%            r2(nl,:) = [X1(nt+1,nl),Y1(nt+1,nl),0];
%            LB(nl) = norm(r2(nl,:));
%         else 
%            r2(nl,:) = [X1(1,nl),Y1(1,nl),0];
%            LB(nl) = norm(r2(nl,:));
%         end
%         CS(nl,:)=cross(r1(nl,:),r2(nl,:));
%         if CS(nl,3)>0
%             sign(nt,nl) = 1;
%         else
%             sign(nt,nl) = -1;
%         end
%         DC(nl)=dot(r1(nl,:),r2(nl,:));
%         if DC(nl)>LA(nl)*LB(nl)
%             DC(nl)=LA(nl)*LB(nl);
%         end
%     end
%     dtheta1(nt,:)=acos(DC./(LA.*LB)).*sign(nt,:);
%     w1(nt,:)=dtheta1(nt,:)/dt;
%     iw1(nt,:)=ic1(nt,:).*w1(nt,:);
%     IW1(nt) = sum(iw1(nt,:));
% end
% for nt=2:num
%     w(nt)=IW1(nt-1)/(Icom1(nt-1));
%     %     rot0(nt)=sum(w(3:nt)*dt);
% end
% 
% for nt = 1:num
%     v_2(nt,:) = vx(nt,:).^2 + vy(nt,:).^2;
%     v(nt,:) = sqrt(v_2(nt,:));
%     vthe(nt,:) = sign(nt,:).*sqrt(v_2(nt,:)-vr(nt,:).^2);  
% end
% for nt = 2:num-1
%     ay(nt,:) = (vy(nt+1,:)-vy(nt-1,:))/(2*dt);
%     ax(nt,:) = (vx(nt+1,:)-vx(nt-1,:))/(2*dt);
%     athe(nt,:) = (vthe(nt+1,:)-vthe(nt-1,:))/(2*dt);
% end
% nt = 1;
% ay(nt,:) = 2*ay(nt+1,:)-ay(nt+2,:);athe(nt,:) = 2*athe(nt+1,:)-athe(nt+2,:);
% ax(nt,:) = 2*ax(nt+1,:)-ax(nt+2,:);
% nt = num;
% ay(nt,:) = 2*ay(nt-1,:)-ay(nt-2,:);athe(nt,:) = 2*athe(nt-1,:)-athe(nt-2,:);
% ax(nt,:) = 2*ax(nt-1,:)-ax(nt-2,:);
% for nt = 1:num
%     c1(nt,:) = 2.*dm.*(X1(nt,:).*vx(nt,:)+Y1(nt,:).*vy(nt,:));
%     c2(nt,:) = dm.*r_2(nt,:);
%     %c3(nt,:) = c1(nt,:).*w(nt,:)+c2(nt,:).*aw(nt,:);
%     %c3(nt,:) = dm.*vr(nt,:).*vthe(nt,:)+dm.*r(nt,:).*athe(nt,:);
%     c3(nt,:) = dm.*(X1(nt,:).*ay(nt,:)-Y1(nt,:).*ax(nt,:));
%     %c3(nt,:) = viw(nt,:);
% end
% % for nt = 1:num
% %     c1(nt,:) = dm.*r_2(nt,:);
% %     c3(nt,:) = dm.*(X1(nt,:).*vy(nt,:)-Y1(nt,:).*vx(nt,:));
% % end
% 
% for nt = 1:num
%     dss(nt,:) = sqrt((x1(nt,2:nl_totl+1)-x1(nt,1:nl_totl)).^2+(y1(nt,2:nl_totl+1)-y1(nt,1:nl_totl)).^2);
% end
% for nt = 1:num
%     C1(nt) = sum(c1(nt,:).*dss(nt,:));
%     C2(nt) = sum(c2(nt,:).*dss(nt,:));
%     C3(nt) = sum(c3(nt,:).*dss(nt,:));
% end
% C3 = C3';
% for nt = 1:num
%     t = nt*dt;
%     for n = 1:nw
%         cosw(nt,n) = cos(n*omega*t);
%         sinw(nt,n) = sin(n*omega*t);
%         ncosw(nt,n) = n*omega*cos(n*omega*t);
%         nsinw(nt,n) = n*omega*sin(n*omega*t);
%     end
% end
% for nt = 1:num
%     t = nt*dt;
%     A(nt,1:nw) = C1(nt)*cosw(nt,:) - C2(nt)*nsinw(nt,:);
%     A(nt,nw+1:2*nw) = C1(nt)*sinw(nt,:) + C2(nt)*ncosw(nt,:);
% %     A(nt,1:nw) = C1(nt)*cosw(nt,:);
% %     A(nt,nw+1:2*nw) = C1(nt)*sinw(nt,:);
% 
% end 
% b0 = [-0.2793 0 0.0016 0 0 0 2.5340 0 -0.0007 0 0 0] ;
% 
% 
% options = optimset('MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-10);
% 
% [coe,result,~,report] = fminsearch(@(B)myfun(A,B,C3,nw,dt),b0,options);
% 
%  for n = 1:num
%     wb(n)=0;
%     for nn= 1:nw
%         wb(n) = wb(n) + coe(nn)*cosw(n,nn);
%     end
%     for nn=1:nw
%         wb(n) = wb(n) + coe(nn+nw)*sinw(n,nn);
%     end
% end
% 
% rot0=cumsum(wb)*dt;
% rot=rot0-mean(rot0);drot=diff(rot)/dt;ddrot=diff(drot)/dt;
% for nt = 1:num
%     x1(nt,:)=x1(nt,:)-cx1(nt);
%     y1(nt,:)=y1(nt,:)-cy1(nt);    
%     x2(nt,:) = x1(nt,:)*cos(rot(nt)) - y1(nt,:)*sin(rot(nt));
%     y2(nt,:) = y1(nt,:)*cos(rot(nt)) + x1(nt,:)*sin(rot(nt));   
%     X2(nt,:) = X1(nt,:)*cos(rot(nt)) - Y1(nt,:)*sin(rot(nt));
%     Y2(nt,:) = Y1(nt,:)*cos(rot(nt)) + X1(nt,:)*sin(rot(nt));    
%     cx2(nt)=sum(dm.*X2(nt,:))/m;
%     cy2(nt)=sum(dm.*Y2(nt,:))/m;
% end
% for nt=1:num
%     X2(nt,:)=X2(nt,:)-cx2(nt);
%     Y2(nt,:)=Y2(nt,:)-cy2(nt);
%     
%     x2(nt,:)=x2(nt,:)-cx2(nt);
%     y2(nt,:)=y2(nt,:)-cy2(nt);
%     
%     ic2(nt,:)=dm.*(X2(nt,:).^2+Y2(nt,:).^2);
%     Icom2(nt)=sum(ic2(nt,:));
% end
for i=2:num-1
%     DI(i)=(Icom2(i+1)-Icom2(i-1))/(2*dt);
     DI(i)=(Icom1(i+1)-Icom1(i-1))/(2*dt);
end
DI(1)=2*DI(2)-DI(3);
DI(num)=2*DI(num-1)-DI(num-2);

num_sp = length(s)-2;
num_larva = 40;
la =  0.02; lb = 0.02;
for m = 1:num_larva
    load('tri0.mat')
    for n = 1:num_sp
        by(n,m) = la*cos((m-1)/num_larva*2*pi);
        bz(n,m) = lb*sin((m-1)/num_larva*2*pi);
        point((n-1)*num_larva+m, 1:3) = [s(n+1), by(n,m), bz(n,m)];
    end
end
np = length(point);
pnt_front(1,1:3) = [s(1) 0 0];
pnt_front(2:np+1,1:3) = point;
pnt_front(np+2,1:3) = [s(end) 0 0];
figure
trimesh(tri(:,:),pnt_front(:,1), pnt_front(:,2),pnt_front(:,3));
xlabel('X','FontName','Times','FontSize',20);
ylabel('Y','FontName','Times','FontSize',20);
zlabel('Z','FontName','Times','FontSize',20);
axis equal
axis tight
% 
for nt = 1:num
    point2(:,3)= point(:,3);
%     figure
%     plot(x1(nt,:),y1(nt,:))
    for i = 1:num_sp
        alpha_tg(nt,i) = atan((Y1(nt,i+1)-Y1(nt,i))/(X1(nt,i+1)-X1(nt,i)));
    end
    i = 1;
%     for i = 1:num_sp
     for j = 1:num_larva
     xtemp=-point((i-1)*num_larva+j,2)*sin(alpha_tg(nt,i))+x1(nt,i+1) ;
     ytemp= point((i-1)*num_larva+j,2)*cos(alpha_tg(nt,i))+y1(nt,i+1);
     point2(num_larva*(i-1)+j,1)=xtemp;
     point2(num_larva*(i-1)+j,2)=ytemp;
%         end
    end
     for i = 2:num_sp
        if ((abs(alpha_tg(nt,i-1))/alpha_tg(nt,i-1))*(abs(alpha_tg(nt,i))/alpha_tg(nt,i))>0 ||((abs(alpha_tg(nt,i-1))/alpha_tg(nt,i-1))*(abs(alpha_tg(nt,i))/alpha_tg(nt,i))<0&&abs(alpha_tg(nt,i))<1.5))            
          for j=1:num_larva
            xtemp=-point((i-1)*num_larva+j,2)*sin(alpha_tg(nt,i))+x1(nt,i+1) ;
            ytemp= point((i-1)*num_larva+j,2)*cos(alpha_tg(nt,i))+y1(nt,i+1);
            point2(num_larva*(i-1)+j,1)=xtemp;
            point2(num_larva*(i-1)+j,2)=ytemp;
          end
        elseif ((abs(alpha_tg(nt,i-1))/alpha_tg(nt,i-1))*(abs(alpha_tg(nt,i))/alpha_tg(nt,i))<0&&abs(alpha_tg(nt,i))>1.5)
            ii1 = i;
            break
        end 
     end
     a1 = exist('ii1');
     if a1 == 1
     if (ii1<num_sp || ii1 == num_sp)
     i = ii1;
     for j = 1:num_larva
     xtemp=point((i-1)*num_larva+j,2)*sin(alpha_tg(nt,i))+x1(nt,i+1) ;
     ytemp=-point((i-1)*num_larva+j,2)*cos(alpha_tg(nt,i))+y1(nt,i+1);
     point2(num_larva*(i-1)+j,1)=xtemp;
     point2(num_larva*(i-1)+j,2)=ytemp;
%         end
     end
    if ii1<num_sp
    for i = ii1+1:num_sp
        if ((abs(alpha_tg(nt,i-1))/alpha_tg(nt,i-1))*(abs(alpha_tg(nt,i))/alpha_tg(nt,i))>0 ||((abs(alpha_tg(nt,i-1))/alpha_tg(nt,i-1))*(abs(alpha_tg(nt,i))/alpha_tg(nt,i))<0&&abs(alpha_tg(nt,i))<1.5))                     
        for j=1:num_larva
            xtemp=point((i-1)*num_larva+j,2)*sin(alpha_tg(nt,i))+x1(nt,i+1) ;
            ytemp=-point((i-1)*num_larva+j,2)*cos(alpha_tg(nt,i))+y1(nt,i+1);
            point2(num_larva*(i-1)+j,1)=xtemp;
            point2(num_larva*(i-1)+j,2)=ytemp;
        end
        elseif ((abs(alpha_tg(nt,i-1))/alpha_tg(nt,i-1))*(abs(alpha_tg(nt,i))/alpha_tg(nt,i))<0&&abs(alpha_tg(nt,i))>1.5)
            ii2 = i;
            break
        end 
    end
    end
    end
    end
    a2 = exist('ii2');
    if a2 == 1
    if (ii2<num_sp || ii2 == num_sp)
    i = ii2;
     for j = 1:num_larva
     xtemp=-point((i-1)*num_larva+j,2)*sin(alpha_tg(nt,i))+x1(nt,i+1) ;
     ytemp=point((i-1)*num_larva+j,2)*cos(alpha_tg(nt,i))+y1(nt,i+1);
     point2(num_larva*(i-1)+j,1)=xtemp;
     point2(num_larva*(i-1)+j,2)=ytemp;
%         end
     end
    if ii2<num_sp
    for i = ii2+1:num_sp
        if ((abs(alpha_tg(nt,i-1))/alpha_tg(nt,i-1))*(abs(alpha_tg(nt,i))/alpha_tg(nt,i))>0 ||((abs(alpha_tg(nt,i-1))/alpha_tg(nt,i-1))*(abs(alpha_tg(nt,i))/alpha_tg(nt,i))<0&&abs(alpha_tg(nt,i))<1.5))                        
        for j=1:num_larva
            xtemp=-point((i-1)*num_larva+j,2)*sin(alpha_tg(nt,i))+x1(nt,i+1) ;
            ytemp=point((i-1)*num_larva+j,2)*cos(alpha_tg(nt,i))+y1(nt,i+1);
            point2(num_larva*(i-1)+j,1)=xtemp;
            point2(num_larva*(i-1)+j,2)=ytemp;
        end
        elseif ((abs(alpha_tg(nt,i-1))/alpha_tg(nt,i-1))*(abs(alpha_tg(nt,i))/alpha_tg(nt,i))<0&&abs(alpha_tg(nt,i))>1.5)
            ii3 = i;
            break
        end 
    end 
    end
    end
    end 
    clear ii1 ii2 ii3
    pnt_front2(nt,1,1:3)=[x1(nt,1),y1(nt,1),0];
    pnt_front2(nt,2:np+1,1:3) = point2;
%     pnt_front2(nt,2:num_larva+1,3) = point2(1:num_larva,3)/3;
%     pnt_front2(nt,num_larva+2:2*num_larva+1,3) = point2(num_larva+1:2*num_larva,3)*2/3;
%     pnt_front2(nt,(num_sp-2)*num_larva+2:(num_sp-1)*num_larva+1,3) = point2((num_sp-2)*num_larva+1:(num_sp-1)*num_larva,3)*2/3;
%     pnt_front2(nt,(num_sp-1)*num_larva+2:np+1,3) = point2((num_sp-1)*num_larva+1:np,3)/3;
    pnt_front2(nt,np+2,1:3)=[x1(nt,100),y1(nt,100),0];    
end
func_inputfiles_creation(pnt_front2);
%  set(0,'defaultfigurecolor','w')
% skip = 10;
%figure;

% for nt = 1:skip: num
%     nt
%     figure
%     trimesh(tri,pnt_front2(nt,:,1), pnt_front2(nt,:,2),pnt_front2(nt,:,3));
% %         trimesh(tri,pnt(nt,:,1), pnt(nt,:,2),pnt2(nt,:,3));
% 
%     %trimesh(tri,ll(nt,:), pnt_f(nt,:,2),pnt_f(nt,:,3), value2);
%     xlabel('X','FontName','Times','FontSize',20);
%     ylabel('Y','FontName','Times','FontSize',20);
%     zlabel('Z','FontName','Times','FontSize',20);
% %     axis([-0.5 0.5 -0.5 0.5 -0.03 0.03])
%     axis equal
%     i=num2str(nt);
%     saveas(gcf,[i,'.png']);
%     %pause
% end
function f = myfun(A,B,C3,nw,dt)
f = 0;
for n = 1:2*nw
    f = f + A(:,n)*B(n);
end
f = f + C3;
f = sum(abs(f))*dt;
end

function [] = plot1( x,y,r )
theta=0:0.1:2*pi;
Circle1=x+r*cos(theta);
Circle2=y+r*sin(theta);
c=[123,14,52];
plot(Circle1,Circle2,'b','linewidth',1);
axis equal
end
