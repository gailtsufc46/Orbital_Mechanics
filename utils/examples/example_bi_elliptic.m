%% Consider the bi-elliptic Hohmann transfer

% The 3 delta-v's for a bi-elliptic maneuver, scaled by sqrt(mu/rA)
dv1f = @(a,b) sqrt( 2*b./(1+b) ) - 1;
dv2f = @(a,b) abs(sqrt( 2*a./b./(a+b) ) - sqrt(2./b./(1+b)));
dv3f = @(a,b) abs(sqrt( 2*b./a./(a+b) ) - sqrt(1./a));

dvtotbef = @(a,b) dv1f(a,b) + dv2f(a,b) + dv3f(a,b);

% the 2 delta-v's for a Hohman transfer, scaled by sqrt(mu/rA)
dvh1f = @(a) sqrt(2*a./(1+a)) - 1;
dvh2f = @(a) sqrt(1./a) - sqrt(2./a./(1+a));

dvtothf = @(a) dvh1f(a) + dvh2f(a);

%% Case 1: rC = 2*rA ... rB = 4*rA
alpha = 2;
beta = 4;

dvtotbef(alpha,beta)
PlotBiElliptic(1,beta,alpha,1)

dvtothf(alpha)
PlotHohmann(1,alpha,1)


%% Case 2: rC = 20*rA ... rB = 40*rA
alpha = 20;
beta = 40;

dvtotbef(alpha,beta)
PlotBiElliptic(1,beta,alpha,1)

dvtothf(alpha)
PlotHohmann(1,alpha,1)



%% Compare total delta-v for Hohmann and Bi-Elliptic

alpha = linspace(2,20);
beta = [10 15 20 40 60];
beta = [10, 11.94, 15.58, 20];
beta = sort([ logspace(0.5,4,100), 11.94, 15.58]);

figure, 
%plot(alpha,dvh1f(alpha),alpha,dvh2f(alpha),alpha,dvtothf(alpha),'linewidth',2)
plot(alpha,dvtothf(alpha),'g','linewidth',3)
grid on, hold on
xlabel('Alpha = r_C / r_A')
ylabel('Total Normalized \DeltaV')
legend('Hohmann')

nalp = 1e4;
dvTotBE = zeros(length(beta),nalp);
for j=1:length(beta)
  alp = linspace(1.5,beta(j),nalp);
  dvTotBE(j,:) = dvtotbef(alp,beta(j));
  if( beta(j) == 15.58 )
    plot(alp,dvTotBE(j,:),'r','linewidth',3,...
      'displayname',sprintf('%s = %1.2f','BE, \beta',beta(j)))
  elseif( beta(j) == 11.94 )
    plot(alp,dvTotBE(j,:),'b','linewidth',3,...
      'displayname',sprintf('%s = %1.2f','BE, \beta',beta(j)))
  elseif( beta(j) < 11.94 )
    plot(alp,dvTotBE(j,:),'c','linewidth',.5,...
      'handlevisibility','off')
  elseif( beta(j) > 15.58 )
    plot(alp,dvTotBE(j,:),'m','linewidth',.5,...
      'handlevisibility','off')
  end
  
end
set(gca,'xlim',[0 alpha(end)])


