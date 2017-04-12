function [P]=gen_mat_P(corner, image_pt)
    [x,y] = size(corner);
    corner = [corner 1];
    Comp_zeros = zeros(1,(y+1));
    P = zeros(2,3*(y+1));
    P(1,:) = [corner Comp_zeros -image_pt(1,1)*corner];
    P(2,:) = [Comp_zeros corner -image_pt(1,2)*corner];
