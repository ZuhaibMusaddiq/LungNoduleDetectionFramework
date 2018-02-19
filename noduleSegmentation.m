function Nodules=noduleSegmentation(lungs,iter)
u = demo_acwe(lungs,iter);
       seg = u<=0;
erodedBW=seg;

   for ukp=0:45:180
        se=strel('line',2,ukp);
   erodedBW = imerode(erodedBW,se);
   end
   
   Nodules=erodedBW;
end