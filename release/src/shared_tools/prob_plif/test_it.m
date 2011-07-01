load plif_data
bins=50 ;

L2 = (L==1 &( B==A | B==Z))*2-1 ;

[prob, prob_cum,limits,prob_unsort,prob_cum_unsort] = ...
    learn_plifs_simple(Out,L2,bins);

[prob, prob_cum,limits,prob_unsort,prob_cum_unsort] = ...
    learn_plifs(Out,L,B,A,Z,bins);

[prob, prob_cum,limits,prob_unsort,prob_cum_unsort] = ...
    learn_plifs_qp(Out,L2,bins);

figure(1) ; clf
plot(prob_unsort) ;
hold on 
plot(prob,'r')

figure(2) ; clf
plot(prob_cum_unsort) ;
hold on 
plot(prob_cum,'r')

figure(3) ; clf
plot(prob) ;
hold on 
plot(prob_cum,'r')


