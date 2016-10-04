close all;

% figure;
% 
% x = -1000:10:1000;
% ph = 1:0.5:14;
% for n=1:length(ph)
%     z = find(phyce.pHRange == ph(n));
%     hist(phyce.charge(z,:),x);
%     xlabel('Charge');
%     ylabel('Counts');
%     title(strcat('pH = ',int2str(ph(n))));
%     M(n) = getframe;
% end
% 
% numtimes = 1;
% fps = 1;
% movie(M,numtimes,fps);

figure(1);
filename = 'phrangeSmallAcidic.gif';
load phyce.mat;

x = -1000:10:1000;
ph = 4:0.2:7;
for n=1:length(ph)
    z = find(phyce.pHRange == ph(n));
    hist(phyce.charge(z,:),x);
    axis([-1000 1000 0 600]);
    xlabel('Charge');
    ylabel('Counts');
    title(strcat('pH = ',num2str(ph(n))));
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if n==1
        imwrite(imind,cm,filename,'gif','Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    %M(n) = getframe;
end

% numtimes = 1;
% fps = 1;
% movie(M,numtimes,fps);