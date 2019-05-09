function [mic_dis, azi_2d, ele_2d] = system_setup(d, res)
    
    x_loc = 0;
    y_loc = 0;
    c = 343;
    ang = [120 60 0 300 240 180];
    mic_loc = [x_loc+d*cosd(ang') y_loc+d*sind(ang')];
    
    mic_dis = [];
    for ref_mic = 1:6
        for idx=1:6
            dis_vec = [mic_loc(idx,1)-mic_loc(ref_mic,1) mic_loc(idx,2)-mic_loc(ref_mic,2) 0];
            mic_dis = [mic_dis; dis_vec];
        end
    end
    
    azi_grid = linspace(0,360,360/res);
    ele_grid = linspace(0,90,90/res);
    
    [azi_2d, ele_2d] = meshgrid(azi_grid,ele_grid);

end