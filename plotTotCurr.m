function plotTotCurr

figure(11)
load( totCurr )

subplot( 2, 1, 1 )
plot( totCurr(:,3)./(totCurr(:,2)+totCurr(:,3)), totCurr(:,1) )
subplot( 2, 1, 2 )
plot( totCurr(:,3)./(totCurr(:,2)+totCurr(:,3)), (totCurr(:,1)-totCurr(1,1))/totCurr(1,1) )
set( gca, 'YLim', 0.01*[-1 1] )
set( gca, 'XLim', [1 size( totCurr, 1 )] )
hold on
plot( totCurr(:,3)./(totCurr(:,2)+totCurr(:,3)), (totCurr(:,1)-totCurr(1,1))/totCurr(1,1), 'r*' )


