function [pointCoord] = findLocation(angles, arrays, d, plotting, sig_name)
    alpha1 = angles(1);
    alpha2 = angles(2);
    alpha3 = angles(3);

    t1=100.00;
    t2=-100.00;
    shift = [-d/2 -sqrt(3)*d/2; d/2 -sqrt(3)*d/2; -d 0; d 0; -d/2 sqrt(3)*d/2; d/2 sqrt(3)*d/2];
    mic1 = repmat([arrays(1,1) arrays(1,2)],6,1); %[238 490; 238 490; 238 490; 238 490; 238 490; 238 490;];
    mic2 = repmat([arrays(2,1) arrays(2,2)],6,1); %[124 338; 124 338; 124 338; 124 338; 124 338; 124 338;];
    mic3 = repmat([arrays(3,1) arrays(3,2)],6,1); %[293 130; 293 130; 293 130; 293 130; 293 130; 293 130;];

    mic1 = mic1+shift;
    mic2 = mic2+shift;
    mic3 = mic3+shift;

    line1point1 = [arrays(1,1) arrays(1,2)]; 
    line2point1 = [arrays(2,1) arrays(2,2)];
    line3point1 = [arrays(3,1) arrays(3,2)];
  
    %%%%% second points for lines
    line1p2 = [line1point1(1)+0.10*cosd(alpha1) line1point1(2)+0.10*sind(alpha1)];
    line2p2 = [line2point1(1)+0.10*cosd(alpha2) line2point1(2)+0.10*sind(alpha2)];
    line3p2 = [line3point1(1)+0.10*cosd(alpha3) line3point1(2)+0.10*sind(alpha3)];

    %%%line Vector
    vector1=line1p2-line1point1;
    vector2=line2p2-line2point1;
    vector3=line3p2-line3point1;

    %two points  for t1 and t2
    %---line 1---
    line1r1 = line1point1 + t1*vector1;
    line1r2 = line1point1 + t2*vector1;

    %---line 2---
    line2r1 = line2point1 + t1*vector2;
    line2r2 = line2point1 + t2*vector2;

    %---line 3---
    line3r1 = line3point1 + t1*vector3;
    line3r2 = line3point1 + t2*vector3;


    %%%%%%%%%%%%%%%%
    %%Find intersection: p for Rp = q
    %%%%%%%%%%%%%%%%
    vector1=(vector1/norm(vector1))';
    vector2=(vector2/norm(vector2))';
    vector3=(vector3/norm(vector3))';
    R = (eye(2)-vector1*vector1')+(eye(2)-vector2*vector2')+(eye(2)-vector3*vector3');
    q =(eye(2)-vector1*vector1')*line1point1'+(eye(2)-vector2*vector2')*line2point1'+(eye(2)-vector3*vector3')*line3point1';

    pointCoord = linsolve(R,q);
    
    %%%%%%%%%%%%%%%%
    %%Plot Lines 
    %%%%%%%%%%%%%%%%
    if(plotting==1)
        fig = figure;
        scatter(mic1(:,1),mic1(:,2), 20, 'fill', 'MarkerFaceColor',[0.6350 0.0780 0.1840]);
        hold on
        scatter(mic2(:,1),mic2(:,2), 20, 'fill', 'MarkerFaceColor', 'b');
        scatter(mic3(:,1),mic3(:,2), 20, 'fill', 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
        scatter(line1point1(1),line1point1(2), 24, 'fill', 'MarkerFaceColor',[0.6350 0.0780 0.1840]);
        scatter(line2point1(1),line2point1(2), 24, 'fill', 'MarkerFaceColor','b');
        scatter(line3point1(1),line3point1(2),24, 'fill', 'MarkerFaceColor',[0.4940 0.1840 0.5560]);

        %%plot lines:
        hold on
        plot([line1r1(1),line1r2(1)], [line1r1(2),line1r2(2)], 'Color', [0.6350 0.0780 0.1840]);

        plot([line2r1(1),line2r2(1)], [line2r1(2),line2r2(2)], 'b');
        plot([line3r1(1),line3r2(1)], [line3r1(2),line3r2(2)],'Color', [0.4940 0.1840 0.5560]);
        %%add point to plot
        str = {'A1','A2','A3',['(' num2str(pointCoord(1)) ', ' num2str(pointCoord(2)) ')']};
        scatter(pointCoord(1), pointCoord(2), 24, 'fill', 'MarkerFaceColor',[0.4660 0.6740 0.1880]);
        %text([arrays(1,1)-5 arrays(2,1)-5 arrays(3,1) pointCoord(1)+5], [arrays(1,2)+15 arrays(2,2)+15 arrays(3,2)+15 pointCoord(2)+2], str);
        xlim([0.00 6.00]);
        ylim([0.00 6.00]);
        title([sig_name '- ' 'User Location: ' '(' num2str(pointCoord(1)) ', ' num2str(pointCoord(2)) ')']);
        hold off
        saveas(fig, ['./results/' sig_name '.png']);
        
    end
end