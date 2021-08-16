% Function Name: gen_net.m

% Authors: David W. Johnson, Tuhin K. Roy and Timothy W. Secomb

% Version History:
%   1.0 - Initial creation

% Description
%   (1) Purpose: Automate resistnet to create three figures: 
%       (A) Cumulative resistance graph for specified ranges with constant
%       constriction (0, 15, 30, 45%) - Figure 3
%       (B) Cumulative resistance graph for specified ranges with variable
%       constriction (0, 15, 30, 45%) - Figure 4A and 4B
%       (C) Cumulative resistance graph for specified ranges with variable
%       constriction (0, 15, 30, 45%) with variable branching ratios -
%       Figure 2
%   (2) Algorithms or Techniques: N/A

% Input
%   None. Inputs are specified in gen_net.m for use in resistnet.m

% Output 
%   Figures (A), (B), (C) as noted above

function [network_resistance] = gen_net()


% (A) Cumulative resistance graph for specified ranges with constant 
% constriction (0, 15, 30, 45%) (AKA Figure 3)
n = 0;
for sqrange1 = 1:6
    deltaresistance = NaN;
    ldratio_correction = NaN;
    nn = 0;
    sqrange2 = sqrange1;
    absvalsqueeze = 5;
    figure
    
    for sqpct1 = 100:-15:55
        sqpct2 = sqpct1;
        nn = 1;
        n = n+1;
        tic
        [~,L3_resistance,ldratio_correction] = resistnet(2, 3, log10(3808),0,sqrange1,sqpct1,sqrange2,sqpct2,nn,1,deltaresistance,absvalsqueeze,ldratio_correction,0);
        network_resistance(:,:,n) = L3_resistance(:,:);
        toc
        set(gcf, 'Position',  [100, 100, 600, 600])

            if sqrange1 == 1
                a = '0 um';
                a = '00 um';
                xlow = 0;
            elseif sqrange1 == 2
                a = '25 um';
                a = '20 um';
                xlow = log10((20^3)/6^3);
            elseif sqrange1 == 3
                a = '50 um';
                a = '60 um';
                xlow = log10((60^3)/6^3);
            elseif sqrange1 == 4
                a = '100 um';
                a = '180 um';
                xlow = log10((180^3)/6^3);
            elseif sqrange1 == 5
                a = '200 um';
                a = '540 um';
                xlow = log10((540^3)/6^3);
            elseif sqrange1 == 6
                a = '400 um';
                a = '1620 um';
                xlow = log10((1620^3)/6^3);
            elseif sqrange1 == 7
                a = '800 um';
                a = '4860 um';
                xlow = log10((4860^3)/6^3);
            elseif sqrange1 == 8
                a = '1 cm';
            end

            if sqrange2 == 1
                b = '25 um';
                b = '20 um';
                xhigh = log10((20^3)/6^3);
            elseif sqrange2 == 2
                b = '50 um';
                b = '60 um';
                xhigh = log10((60^3)/6^3);
            elseif sqrange2 == 3
                b = '100 um';
                b = '180 um';
                xhigh = log10((180^3)/6^3);
            elseif sqrange2 == 4
                b = '200 um';
                b = '540 um';
                xhigh = log10((540^3)/6^3);
            elseif sqrange2 == 5
                b = '400 um';
                b = '1620 um';
                xhigh = log10((1620^3)/6^3);
            elseif sqrange2 == 6
                b = '800 mm';
                b = '4860 um';
                xhigh = log10((4860^3)/6^3);
            elseif sqrange2 == 7
                b = '1 cm';
            elseif sqrange2 == 8
                b = '10 cm';
            end
            
    end
    
    txt = ['A'];
    text(1,180,txt,'FontSize',30, 'HorizontalAlignment','left','FontWeight','bold')

    figuretitle = char(strcat('Fig.3_[',a,'_',b,']','_[',string(sqpct1),'%_',string(sqpct2),'%].fig'));
    savefig(figuretitle);

    %Create shaded boxes
    hold on
    v = [xlow 0; xhigh 0; xhigh 200; xlow 200];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices',v,'FaceAlpha', 0.025,'EdgeColor','none');
    hold off
    
end

clear all

%Fig 4A -  HPV Best Guess Simulation with upward ramp (log(cap) dependent),
%plateau from 100-400um and downward ramp and rb = 3
nn = 0;
ldratio_correction = [NaN];
sqrange1 = [NaN];
sqrange2 = [NaN];
deltaresistance = [NaN];
absvalsqueeze = 2;                      %2 is best guess setting
figure
for sqpct1 = 100:-15:55     %% Put in variable for loop
    sqpct2 = sqpct1;
    nn = 1;
    tic
    [~,L3_resistance,ldratio_correction] = resistnet(2, 3, log10(3808),1,sqrange1,sqpct1,sqrange2,sqpct2,nn,1,deltaresistance,absvalsqueeze,ldratio_correction,0);
    network_resistance(:,:,nn) = L3_resistance(:,:);
    toc
        
        set(gcf, 'Position',  [100, 100, 600, 600])
        ylim([0 200])

        if sqpct1 == 100 
            %Mild resistance elevation from hypoxia (less effects of increased CO) in humans (baseline multiplied by factor)
            %Create HPV Reference from 100% diam from Talbot
            %https://journals.physiology.org/doi/pdf/10.1152/japplphysiol.00903.2004
            if sqpct1 == 100
                hpvx = [10.74];
                %56 is est HPV elevation see below
                raw = [1.8];
                factor = raw./1.1;
                
                hpvresistance = L3_resistance(3,6)*factor;
                hpvy = [hpvresistance];
                hold on
                plot(hpvx,hpvy,'kp','MarkerSize',15)
                hold off
            
            
            %Spectrum of resistance elevation factor in canines from zero to complete hypoxia
            %Create HPV Reference from 100% Diam from Marshall 1994
            %https://www-sciencedirect-com.ezproxy1.library.arizona.edu/science/article/pii/0034568794901295
                hpvx = [10.74 10.74 10.74];
                %Order of raw is 100mmHg O2, 50, 0
                raw = [1.2 1.8 3.3];
                factor = raw./1.2;
                letterlabel = {'A' 'B' 'C'};
                
                %56 is est HPV elevation see below
                hpvresistance = L3_resistance(3,6)*factor;
                hpvy = [hpvresistance];
                hold on
                for i = 1:3
                    text(hpvx(i), hpvy(i), letterlabel(i), 'HorizontalAlignment','center', 'VerticalAlignment','middle','FontSize',15,'FontWeight','bold')
                end
                %plot(hpvx,hpvy,'kd')
                hold off
            
%             %Create Elevation of PVR with exercise (same as normal)
%             %https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3204788/
%             %PVR remains constant at 1.2 Wood units (mmHg/L/min) w/ exercise
%             %1.2 Wood units = 1.2 (80) Dynes sec cm^-5 = 96
%                 hpvx = [10.74];
%                 %56 is est HPV elevation see below
%                 hpvresistance = 96;
%                 hpvy = [hpvresistance];
%                 hold on
%                 plot(hpvx,hpvy,'gs')
%                 hold off
            end 
            
        end
        savefig('Fig.4_Best_Guess.fig');   
end
nn = 0;


% Figure 4A - Profile of vasoconstriction for corresponding to 4B
figure
for sqpct = 0:15:45
    
    xguess = [0 33  100    300      900  1000];
    yguess = [0 0   sqpct  sqpct    0       0];
    
    totallines = 100;
    colmap = hsv(totallines);

    hold on
    plot(xguess,yguess,'Color',colmap((100-sqpct),:))
    hold off
    
    txt = [num2str(sqpct) '%'];
    text(180,sqpct+1,txt,'FontSize',15)
    
    xlabel('Vessel Diameter (um)', 'FontSize', 20)
    ylabel('Percent Constriction (%)', 'FontSize', 20)
    set(gcf, 'Position',  [100, 100, 600, 600])
    ylim([0 1000])
    ylim([0 50])
    
    savefig('Fig.4_Best_Guess_PctProfile.fig');
    
end

clear all

%Fig 2 - Example of the effects of vasoconstriction at multiple ranges of 
% vasoconstriction (30% constriction HPV simulation over one diameter 
% range with varied branching ratio.
for sqrange1 = 1:6
    deltaresistance = NaN;
    ldratio_correction = NaN;
    nn = 0;
    sqrange2 = sqrange1;
    absvalsqueeze = 0;
    
    for sqpct1 = 100:-30:70
        sqpct2 = sqpct1;
        nn = 1;
        tic
        [~,L3_resistance,ldratio_correction] = resistnet(2, 4, log10(3808),0,sqrange1,sqpct1,sqrange2,sqpct2,nn,1,deltaresistance,absvalsqueeze,ldratio_correction,0);
        network_resistance(:,:,nn) = L3_resistance(:,:);
        toc
        set(gcf, 'Position',  [100, 100, 600, 600])

            if sqrange1 == 1
                a = '0 um';
                a = '20 um';
            elseif sqrange1 == 2
                a = '25 um';
                a = '20 um';
            elseif sqrange1 == 3
                a = '50 um';
                a = '60 um';
            elseif sqrange1 == 4
                a = '100 um';
                a = '180 um';
            elseif sqrange1 == 5
                a = '200 um';
                a = '540 um';
            elseif sqrange1 == 6
                a = '400 um';
                a = '1620 um';
            elseif sqrange1 == 7
                a = '800 um';
                a = '4860 um';
            elseif sqrange1 == 8
                a = '1 cm';
            end

            if sqrange2 == 1
                b = '25 um';
                b = '20 um';
            elseif sqrange2 == 2
                b = '50 um';
                b = '60 um';
            elseif sqrange2 == 3
                b = '100 um';
                b = '180 um';
            elseif sqrange2 == 4
                b = '200 um';
                b = '540 um';
            elseif sqrange2 == 5
                b = '400 um';
                b = '1620 um';
            elseif sqrange2 == 6
                b = '800 mm';
                b = '4860 um';
            elseif sqrange2 == 7
                b = '1 cm';
            elseif sqrange2 == 8
                b = '10 cm';
            end

            figuretitle = char(strcat('Fig.5_[',a,'_',b,']','_[',string(sqpct1),'%_',string(sqpct2),'%].fig'));

            %Extract data from a graph
            %x=get(h,'Xdata')
            %y=get(h,'Ydata')

    end
    savefig(figuretitle);
    nn = 0;
end


clear all

%Finish notification
sound(sin(1:3000))
sound(sin(1:3000))
sound(sin(1:3000))
end