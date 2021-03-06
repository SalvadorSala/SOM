
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load modelexample

% close all
fila=6 % column in the model that cuantifies the number fo correct maps for each iteration
hasta=133;
ejes=[0 11 50 100]
Tamletra=20;
ticks=[0,0];

h2=figure('Name','','Units','default')
h2Title  = title ('');
h2XLabel = xlabel('\sigma_{molecular}');
h2YLabel = ylabel('%Correct Orientation');
h2Legend = legend('0.01','0.03','0.06');
ax1=gca;
axes(ax1)
inicio=1:44:hasta-1;
sig_ephrin=[0.1,1:10]
alfan=0:0.03:0.21;
counter=0;
colr1={[0.1333    0.5451    0.1333];  [0.6039    0.8039    0.1961]; [0.6784    1.0000    0.1843]};
colr2={ [0         0    0.5451];  [0 0 1]; [ 0.2549    0.4118    0.8824]};
colr3={'red';  [ 0.9100 0.4100 0.1700]; [0.9290 0.6940 0.1250]};

for i=inicio
counter=counter+1;
  ret_wave=[model{i+1:4:i+44,fila}];
hold on
h2plot=plot(sig_ephrin,ret_wave, 'LineWidth',1.5,'color',colr1{counter})
axis(ejes)
h2Legend = legend('0.01','0.03','0.06');
end

inicio=1:44:hasta-1;
sig_ephrin=[0.1,1:10]
alfan=0:0.03:0.21;
counter=0;
for i=inicio
counter=counter+1;
sincro_ini=[model{i+2:4:i+44,fila}];
hold on
h3plot=plot(sig_ephrin,sincro_ini, 'LineWidth',1.5,'color',colr2{counter})
axis(ejes)
h3Legend = legend('0.01','0.03','0.06');
end


inicio=1:44:hasta-1;
sig_ephrin=[0.1,1:10]
alfan=0:0.03:0.21;
counter=0;

for i=inicio
counter=counter+1;
azar=[model{i+4:4:i+44,fila}];
hold on
h4plot=plot(sig_ephrin,azar, 'LineWidth',1.5,'color',colr3{counter})
axis(ejes)
h4Legend = legend('0.01','0.03','0.06');
end



set([ax1], ...
    'FontName'   , 'Arial','LineWidth',1);
set([h2XLabel], ...
    'FontName', 'Arial', 'FontSize', 10,'FontWeight','normal');
set([h2YLabel], ...
    'FontName', 'Arial', 'FontSize', 10,'FontWeight','b');
set([h4Legend], ...
    'FontSize',8,'location', 'SouthWest','FontWeight','bold',...
     'FontName'   , 'Arial');
set([h2Title],...
    'FontSize', 10 ,'FontWeight' , 'bold','FontName', 'Arial','Position',[5.5000 104 0]);
set(gcf,'units','centimeters','position',[10,10, 5.8000 5.8]);
set(gcf,'Paperunits','centimeters','PaperPosition',[1 18  5.8000 5.8])
set([h4Legend],'Visible','off');



set([ax1],'TickDir','out')
set(gcf,'units','normalized')
text(-2,106,'A','FontSize',11,'FontName','Arial','FontWeight','Bold');
set(gcf,'units','centimeters','position',[10,10, 5.8000 5.8]);
set(gcf,'Paperunits','centimeters','PaperPosition',[1 23  5.8000 5.8])

print('NEWFigure3300dpi','-dpdf','-r300')

% fig.PaperPositionMode = 'auto';

set(gcf,'Paperunits','centimeters','PaperOrientation','landscape','PaperPosition',[2 10 26 26])

print('NEWFigure3fillPage300dpi','-dpdf','-r300')
