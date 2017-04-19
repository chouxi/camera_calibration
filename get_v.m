%% Get V matrix
function V=get_v(H)
H=H';
V=[get_v_vec(1,2,H); (get_v_vec(1,1,H)-get_v_vec(2,2,H))];
end

function v=get_v_vec(i,j,H)
v=[H(i,1)*H(j,1) H(i,1)*H(j,2)+H(i,2)*H(j,1) H(i,2)*H(j,2) H(i,3)*H(j,1)+H(i,1)*H(j,3) H(i,3)*H(j,2)+H(i,2)*H(j,3) H(i,3)*H(j,3)];
end
