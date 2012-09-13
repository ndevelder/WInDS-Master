function [sz,sx,vm]=wakeplot(pos,vel,turbine,blade,wake,fastout,j)

opengl autoselect

nb=turbine.NumBl;
ns=size(pos.trail,1);
c=zeros([2*ns+1 3 nb]);

D=2*blade.TipRad;

cc=[0.5 0.5 1 ; 1 0 0 ; 0 1 0];

sx=NaN*ones(ns,1);
sy=sx;
sz=sx;

vx=NaN*ones(ns,1);
vy=vx;
vz=vx;

for k=1:nb
    if k>1
        sx=[sx NaN*ones(ns,1) squeeze(wake.domain{j}(:,1,2:end,k))]; %#ok<AGROW>
        sy=[sy NaN*ones(ns,1) squeeze(wake.domain{j}(:,2,2:end,k))]; %#ok<AGROW>
        sz=[sz NaN*ones(ns,1) squeeze(wake.domain{j}(:,3,2:end,k))]; %#ok<AGROW>
        vx=[vx NaN*ones(ns,1) squeeze(vel.uind{j}(:,1,2:end,k))]; %#ok<AGROW>
        vy=[vy NaN*ones(ns,1) squeeze(vel.uind{j}(:,2,2:end,k))]; %#ok<AGROW>
        vz=[vz NaN*ones(ns,1) squeeze(vel.uind{j}(:,3,2:end,k))]; %#ok<AGROW>
        ct1=[ct1 0*ones(ns,1) cc(k,1)*ones([ns j])]; %#ok<AGROW>
        ct2=[ct2 0*ones(ns,1) cc(k,2)*ones([ns j])]; %#ok<AGROW>
        ct3=[ct3 0*ones(ns,1) cc(k,3)*ones([ns j])]; %#ok<AGROW>
    else
        sx=squeeze(wake.domain{j}(:,1,2:end,k));
        sy=squeeze(wake.domain{j}(:,2,2:end,k));
        sz=squeeze(wake.domain{j}(:,3,2:end,k));
        vx=squeeze(vel.uind{j}(:,1,2:end,k));
        vy=squeeze(vel.uind{j}(:,2,2:end,k));
        vz=squeeze(vel.uind{j}(:,3,2:end,k));
        ct1=cc(k,1)*ones([ns j]);
        ct2=cc(k,2)*ones([ns j]);
        ct3=cc(k,3)*ones([ns j]);
    end
end

vind=mean(sqrt(fastout.WindVxi.^2+fastout.WindVyi.^2+fastout.WindVzi.^2));
vm=sign(vx).*sqrt(vx.^2+vy.^2+vz.^2)/vind;

cb=['g','r','b'];

hold on
for k=1:nb
    s(:,1)=squeeze(wake.domain{j}(end,1,2:end,k));
    s(:,2)=squeeze(wake.domain{j}(end,2,2:end,k));
    s(:,3)=squeeze(wake.domain{j}(end,3,2:end,k));
    plot3(s(:,1)/D,s(:,2)/D,s(:,3)/D,['-' cb(k)],'linewidth',0.5)
end
caxis([-0.1 0.2])
if size(sx,2)>1
    surface(sx/D,sy/D,sz/D,vm,'linewidth',0.025,'linestyle','-')
end
alpha(0.45)
for k=1:nb
    a=squeeze(pos.lead(:,:,j,k));
    b=flipdim(squeeze(pos.trail(:,:,j,k)),1);
    c(:,:,k)=[a ; b ; a(1,:)];    
    fill3(squeeze(c(:,1,k))/D,squeeze(c(:,2,k))/D,squeeze(c(:,3,k))/D,[0.8 0.8 0.8])
    plot3(squeeze(c(:,1,k))/D,squeeze(c(:,2,k))/D,squeeze(c(:,3,k))/D,'-k','linewidth',0.5)
    clear a b s w p v ct
end
fs=10;
set(gca,'fontsize',fs,'fontname','Times New Roman')
hold off
axis tight
axis equal
view([-45 22])
xlabel('Downwind [D]','fontsize',fs,'fontname','Times New Roman')
ylabel('Lateral [D]','fontsize',fs,'fontname','Times New Roman')
zlabel('Vertical [D]','fontsize',fs,'fontname','Times New Roman')
orient landscape