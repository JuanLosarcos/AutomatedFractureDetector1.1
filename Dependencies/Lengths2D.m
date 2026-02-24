function len = Lengths2D(p)
% Lengths2D Computes the Euclidean distance between two 2D points.
% Input: p = [x1 y1 x2 y2]
% Output: scalar length between (x1,y1) and (x2,y2)

    dx = p(3) - p(1);
    dy = p(4) - p(2);
    len = sqrt(dx^2 + dy^2);
end
