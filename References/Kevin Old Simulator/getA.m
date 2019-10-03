function A = getA(ax)
A = zeros(3,3);
for i = 1:length(ax)
    if ax(i) == 'x'
        A(1,i) = 1;
    elseif ax(i) == 'y'
        A(2,i) = 1;
    elseif ax(i) == 'z'
        A(3,i) = 1;
    end
end

end