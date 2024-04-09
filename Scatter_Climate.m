figure
hold on

for i=1:1233
    d=d_1(i,:)*x;
    scatter(d(1),d(2),3,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[0.8 0 0],'LineWidth',0.2);
    if(mod(i,300)==0)
        fprintf("%6.2f percent completed...\n",i*1.0/12332*100);
    end
end
for i=1:668
    d=d_2(i,:)*x;
    scatter(d(1),d(2),3,'MarkerEdgeColor',[0 1 0],'MarkerFaceColor',[0.8 0 0],'LineWidth',0.2);
    if(mod(i,300)==0)
        fprintf("%6.2f percent completed...\n",i*1.0/6688*100);
    end
end