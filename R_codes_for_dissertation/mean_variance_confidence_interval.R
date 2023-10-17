
var(w[,1])

var(w[,2])

var(theta[,1])

var(theta[,2])









paste0("(",round(t.test(w[,1])$conf.int[1],5),  ","  ,round(t.test(w[,1])$conf.int[2],5),")" )
paste0("(",round(t.test(w[,2])$conf.int[1],5),  ","  ,round(t.test(w[,2])$conf.int[2],5),")" )
paste0("(",round(t.test(theta[,1])$conf.int[1],5),  ","  ,round(t.test(theta[,1])$conf.int[2],5),")" )
paste0("(",round(t.test(theta[,2])$conf.int[1],5),  ","  ,round(t.test(theta[,2])$conf.int[2],5),")" )
