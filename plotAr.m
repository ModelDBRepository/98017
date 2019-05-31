% This script plots arrows between the populations in showConnMulti.m

% x_po, y_po: Coordinates of the post-population
% x_pr, y_pr: Coordinates of the pre-population
% an           : Angle of inclination
% x0, y0 : Coordinates of arrow shaft
% x1, y1 : Coordinates of arrow point
% x2, y2 : Coordinates of arrowhead right corner
% x3, y3 : Coordinates of arrowhead left corner

% For auto-arrow
% x4, y4, x5, y5, x6, y6, etc: Other arrow points

%Arrowhead
w2 = pi/6;
r2 = 0.05;
w3 = -pi/6;
r3 = 0.05;

for iiiii = 1:length( arrows )
    x_po = pop( arrows(iiiii).post ).x;  
    y_po = pop( arrows(iiiii).post ).y;
  
    if arrows(iiiii).pre == arrows(iiiii).post % self connection
        d = 0.025;
        an = abs( arrows(iiiii).labpos );
        r = pop( arrows(iiiii).post ).radius;
        [x1,y1] = pol2cart(an, r );
        x7 = x_po + x1 * sqrt( r*r-d*d ) / r;
        y7 = y_po + y1 * sqrt( r*r-d*d ) / r;
        x6 = x7 + 0.1*cos(an);
        y6 = y7 + 0.1*sin(an);
        x1 = x7 + d * cos(an+pi/2);
        y1 = y7 + d * sin(an+pi/2);
        x0 = x1 + 0.1*cos(an);
        y0 = y1 + 0.1*sin(an);
        x4 = x7 - d * cos(an+pi/2);
        y4 = y7 - d * sin(an+pi/2);
        x5 = x4 + 0.1*cos(an);
        y5 = y4 + 0.1*sin(an);
    
       line( [x0 x5], [y0 y5] )
       line( [x4 x5], [y4 y5] )
    
   elseif arrows(iiiii).pre < 0  % external connection
       an = abs( arrows(iiiii).pre );
       x_pr = x_po + 0.1*cos( an );
       y_pr = y_po + 0.1*sin( an );
       [x1,y1] = pol2cart(an, pop( arrows(iiiii).post ).radius );
       x1 = x1 + x_po;
       y1 = y1 + y_po;
       x0 = x1 + 0.1*cos(an);
       y0 = y1 + 0.1*sin(an);
       
       
       
   else  % arrow between populations
       x_pr = pop( arrows(iiiii).pre).x;
       y_pr = pop( arrows(iiiii).pre).y;
       if x_pr == x_po
           an = pi/2+pi*(y_po>y_pr);
       else
           an = atan( (y_pr-y_po)/(x_pr-x_po) );
       end
       if x_pr<x_po
           an = an+pi;
       end
       [x0,y0] = pol2cart(an+pi, pop( arrows(iiiii).pre ).radius );
       x0 = x0 + x_pr;
       y0 = y0 + y_pr;
       [x1,y1] = pol2cart(an, pop( arrows(iiiii).post ).radius );
       x1 = x1 + x_po;
       y1 = y1 + y_po;
   end 
   
    
   
   [x2,y2] = pol2cart(w2+an,r2);
   x2 = x2 + x1;  
   y2 = y2 + y1;
   [x3,y3] = pol2cart(w3+an,r3);
   x3 = x3 + x1;  
   y3 = y3 + y1;
  
   line( [x1 x2], [y1 y2] )
   line( [x1 x3], [y1 y3] )
   line( [x0 x1], [y0 y1] )
    

  
  
   % Plot labels
   txxt(iiiii) =text( 0, 0, arrows(iiiii).label );
   pos = get( txxt(iiiii), 'Extent' );
   wt = pos(3);
   ht = pos(4);
   delete( txxt(iiiii) )
  
      
           
   % d = 1 on right side of arrow, d = -1 on left side of arrow (seen in the
   % direction of the arrow
   d = 2 * arrows(iiiii).labpos - 1;
   
   % kant = 0,1,2,3 when right,upper,left and lower border of text box is next to arrow 
   kant = mod( floor( an/(pi/2) ) + 1 + arrows(iiiii).labpos, 4 );
  
   % dlr = 1 towards left, -1 towards right, dud = 1 upward, -1 downward
   dlr = ( 2 * ( mod( an, 2*pi ) < pi ) - 1 ) * d;
   dud = ( 2 * ( mod( an+pi/2, 2*pi ) < pi ) - 1 ) * d;
   if arrows(iiiii).pre == arrows(iiiii).post % self connection
       xxx(iiiii) = x6;
       yyy(iiiii) = y6; 
       dlr = ( 2 * ( mod( an-pi/2, 2*pi ) < pi ) - 1 );
       dud = ( 2 * ( mod( an, 2*pi ) < pi ) - 1 );
   else
       if d == 1
           xxx(iiiii) = x2;
           yyy(iiiii) = y2;
       else
           xxx(iiiii) = x3;
           yyy(iiiii) = y3;
       end
       if an == 3*pi/2 | an == pi/2 % to get right sign on tan
           an = 3*pi/2+0.000001;
       end
       if mod( kant, 2 )
           if abs( tan(an) ) < 10000
               yyy(iiiii) = y1 + (xxx(iiiii)-x1)*(tan(an));
           else
               yyy(iiiii) = y1 + (xxx(iiiii)-x1)*sign( tan( an ) ) * 10000;
           end
       else
           if abs( tan(an) ) > 1/10000
               xxx(iiiii) = x1 + (yyy(iiiii)-y1)/(tan(an));
           else
               xxx(iiiii) = x1 + (yyy(iiiii)-y1)*1e10;
           end
       end
       if mod( kant+d, 4 ) == 0
           if yyy(iiiii) > y2+ht
               yyy(iiiii) = y2;
               dud = dud*(-1);
           end
       elseif mod( kant+d, 4 ) == 1
           if xxx(iiiii) < x2-wt
               xxx(iiiii) = x2;
               dlr = dlr * (-1);
           end
       elseif mod( kant+d, 4 ) == 2
           if yyy(iiiii) < y2-ht
               yyy(iiiii) = y2;
               dud = dud * (-1);
           end
       elseif mod( kant+d, 4 ) == 3
           if xxx(iiiii) > x2+wt
               xxx(iiiii) = x2;
               dlr = dlr * (-1);
           end
       end
   end
   if dlr > 0
       hstr = 'right';
   else
       hstr = 'left';
   end
   if dud > 0 
       vstr = 'bottom';
   else 
       vstr = 'top';
   end
end

arrrr = arrows;
i = 1;
while length(xxx)>0
    curr = find(xxx == xxx(1) & yyy == yyy(1));
    currlab = [arrrr(curr(1)).label];
    for i = 2:length(curr)
        currlab = [currlab ', ' arrrr(curr(i)).label];
    end
    txxt(i) = text( xxx(1), yyy(1), currlab, ...
 		      'HorizontalAlignment', hstr, 'VerticalAlignment', vstr );
    xxx(curr) = [];
    yyy(curr) = [];
    arrrr(curr) = [];
    i = i + 1;
end
set( gca, 'FontSize', 14 )

   