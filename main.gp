/* x is a generator of Fq obtained via ffgen(q) */
[q,t,m] = readvec("input.txt");
f = ffgen(q,'a); \\création de f, générateur de F128
encodefq(i,x)=subst(Pol(digits(i,x.p)),'x,x);
decodefq(z,x)=fromdigits(Vec(z.pol),x.p);
int2fqx(i,x)=Polrev([encodefq(z,x)|z<-digits(i,x.p^x.f)]);
fqx2int(p,x)=fromdigits([decodefq(z,x)|z<-Vecrev(p)],x.p^x.f);


\\Il s'agit ici simplement d'implémenter les approximants de Padé pour découvrir les erreurs du chiffré, comme pour tout code BCH. On ne sait pas, au départ, à quelle puissance de alpha commencent les racines consécutives, pour cela on introduit une constante de décalage. Enfin, la fonction ffprimroot ne donne pas toujours la meme racine, il faut donc la répéter suffisament souvent pour tomber sur la bonne. On construit alors une liste des erreurs, et si elle est suffisament longue pour pouvoir constituer un message on l'affiche.


y=int2fqx(m, f); 


syndrome(y, a, decal) = sum(i=0, 2 * t - 1, subst(y, 'x, a^(i+decal)) * 'x^i);
bchdecode() = {
	for(k=0, 128,
		a = ffprimroot(f); 
		for(j=0, 127,
			padeapprox=bestapprPade(Mod(syndrome(y, a, j), x^(2*t)));
			Rpade = numerator(padeapprox);
			Epade = denominator(padeapprox); 
			liste_des_erreurs = List(); 
			for(i = 0, 126, pol = subst(Epade, 'x, a^(-i)); if (pol == 0, v = subst(Rpade/deriv(Epade) * (x^(j-1)), 'x, a^(-i)); l = fqx2int(v, f); listput(liste_des_erreurs, l)));
			if(#liste_des_erreurs > 5 ,print(Strchr(Vec(liste_des_erreurs))); return);
		);
	);
};
bchdecode();
