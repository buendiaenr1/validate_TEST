extern crate statrs;

use plotters::style::full_palette::ORANGE;
use plotters::prelude::*;
use plotters::drawing::IntoDrawingArea;
//use plotters::coord::Shift;
use plotters::style::Color;

use statrs::distribution::ChiSquared;
use statrs::assert_almost_eq;
use crate::statrs::distribution::ContinuousCDF;
use statrs::distribution::FisherSnedecor;


use colored::Colorize; // texto de color en msdos

use std::cmp::Ordering;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::{f64, usize};

// ******************************************************************
// autor: enrique ricardo pablo buendia lozada
// 2023 BUAP, Puebla, México
// ******************************************************************

fn f_quantile(p: f64, df1: u64, df2: u64) -> f64 {
    // crear la distribucion FisherSnedecor con df1 y df2
    let f_dist = FisherSnedecor::new(df1 as f64, df2 as f64).unwrap();
    // Usar el método inverse_cdf para calcular el cuantil
    f_dist.inverse_cdf(p)
}
fn suma(vv: Vec<f64>) -> f64 {
    vv.iter().sum()
}
fn elevar_y_sumar(vector: Vec<f64>, potencia: f64) -> f64 {
    let mut suma = 0.0;
    for elemento in vector {
        suma += elemento.powf(potencia);
    }
    suma
}

fn l_linn(x:&Vec<f64>, y:&Vec<f64>) -> f64 {
    // muestras del mismo tamaño
    assert_eq!(x.len(),y.len());
    let n= x.len() as f64;
    // promedios
    let x_m=x.iter().sum::<f64>()/n;
    let y_m=y.iter().sum::<f64>()/n;
    // varianzas
    let x_v=x.iter().map(|&xi|(xi-x_m).powi(2)).sum::<f64>()/n;
    let y_v=y.iter().map(|&yi|(yi-y_m).powi(2)).sum::<f64>()/n;
    // cov
    let cov=x.iter().zip(y.iter()).map(|(&xi, &yi)|(xi - x_m) * (yi - y_m)).sum::<f64>()/n;
    // Pearson
    let cor= cov/(x_v * y_v).sqrt();
    // Lawrence - linn
    let ccc = 2.0 * cor * x_v * y_v / (x_v.powi(2)+ y_v.powi(2) + (x_m - y_m).powi(2));

    ccc
}

fn posthoc_scheffe(medias:&Vec<f64>,cmd:f64,fcritica:f64){
    let mut resultad: Vec<f64> = Vec::new();
    let mut rr: &str;
    let n: usize=medias.len();
    let mut numerosc:Vec<usize>=Vec::new();
    for i in 0..n{
        numerosc.push(i);
    }

    let val_scheffe: f64=((n as f64 -1.0) * fcritica).powf(0.5) * 
        ((1.0/n as f64 + 1.0/n as f64)*cmd).powf(0.5);

    let resp1: &str="_Test Post Hoc Scheffé";
    println!("\n\n{}",resp1.yellow());
    println!(" | {:^15} | {:^15} | {:^15} | {:^15}","vs","Diferencias","Valor Scheffé","Decisión");
    println!(" | {:^15} | {:^15} | {:^15} | {:^15}","--","--","--","--");
    for k in 0..=medias.len()-1 {
        let a = numerosc[k];
        for i in k+1..=medias.len()-1 {
            let b = numerosc[i];
            let cad=a.to_string()+&" ".to_string() +&b.to_string();
            resultad.push((medias[k]-medias[i]).abs());
            let dif: f64=(medias[k]-medias[i]).abs();
            if val_scheffe < dif {
                rr="Significativa";
            }else {
                rr="No Significativa";
            }
            println!(" | {:^15} | {:^15.9} | {:^15.9} | {:^15}",cad,dif,val_scheffe,rr);
        }
    }
}

fn blanalt(x:&Vec<f64>, y:&Vec<f64>,cad:String){
    let resp1 = " Bland -Altman ";
    println!("\n {} ",resp1.yellow());
    assert_eq!(x.len(),y.len());
    //mn<-      (x+y)/2
    let mut mn:Vec<f64> = Vec::new();
    for i in 0..x.len() {
        mn.push((x[i]+y[i])/2.0);
    }
    println!(" >>>>>> promedios {:?}",mn);
    //differ<-  x-y
    let mut differ:Vec<f64> = Vec::new();
    for i in 0..x.len() {
        differ.push(x[i]-y[i]);
    }
    println!(" >>>>>> diferencias {:?}",differ);
    //md<-      mean(differ)
    let n=x.len();
    let mut md=0.0;
    for i in 0..n{
    md= md + differ[i];
    }
    md=md/(n as f64);
    println!(" media de las diferencias {}",md);

    //vd<-      var(differ)
    let mut vd=0.0;
    for i in 0..n{
    vd = vd + (differ[i] - md).powf(2.0);
    }
    vd = vd / (n as f64);

    //plot(mn,differ,pch="*",xlab="media",ylab="diferencia",sub=paste("Media= ",md))
    
    //abline(h=md,lty=2)
    //abline(h=md+c(2,-2)*sqrt(vd),col="red")
    //cat("Intervalo de confidencia del 95%:",md+2*sqrt(vd),"y",md-2*sqrt(vd),"\n")

    let ls = md + 2.0 * vd.sqrt();
    let li = md - 2.0 * vd.sqrt();

    println!(" IC 95% ({},{})",li,ls);

    // Presición de acuerdo a:      
    //   Journal of Clinical Monitoring and Computing
    //           (2008) 22:257–259
    // A MEASURE OF CONFIDENCE
    // IN BLAND–ALTMAN ANALYSIS FOR THE
    // INTERCHANGEABILITY OF TWO METHODS
    // OF MEASUREMENT            
    //        David Preiss, PhD MASc and Joseph Fisher, MD FRCP(C)
    let nn= x.len() as f64;
    // promedios
    let mut x_m=0.0;
    for z in 0..n{
        x_m = x_m + x[z];
    }
    x_m=x_m/nn;
    
    let mut y_m=0.0;
    for z in 0..n{
        y_m = y_m + y[z];
    }
    y_m = y_m/nn;
    
    // varianzas
    let x_v=x.iter().map(|&xi|(xi-x_m).powi(2)).sum::<f64>()/nn;
    let y_v=y.iter().map(|&yi|(yi-y_m).powi(2)).sum::<f64>()/nn;
    let pres = 2.0_f64.sqrt() * (x_v.sqrt() + y_v.sqrt())/2.0;
    let resp1=" precisión ";
    println!(" {} {} mas cercana a cero es mas preciso",resp1.yellow(),pres);
    
    graficar2(&mn, &differ, md, li, ls, cad);
   

        
}

fn m_para_cinco(vector:Vec<f64>) -> usize {

  let mut max:usize=0;  

  for i in 0..vector.len(){
    let cadena: String = vector[i].to_string();
    if cadena.contains(".") { 
        let parts: Vec<&str> = cadena.split ('.').collect::<Vec<&str>> (); // dividir la cadena por el punto decimal
        let decimals: usize = parts [1].len (); // obtener la longitud de la segunda parte
        //println! ("{} tiene {} decimales", cadena, decimals); // imprimir el resultado
        if decimals >= max {
            max=decimals as usize;
        }
    }
  }
  
  return max;
}

fn chi_ccdf(x: f64) -> f64 {
    // los grados de libertad son g.l.=2
    let c: ChiSquared = ChiSquared::new(2.0).unwrap();
    let p: f64 = 1.0 - c.cdf(x);
    assert_almost_eq!(p, 1.0 - c.cdf(x), 1e-15); // verificar si CCDF es igual a 1 - CDF
    p
}


// Define a function to draw a scatter plot with an average line
fn graficar2(xx: &Vec<f64>, yy: &Vec<f64>, med: f64, li: f64, ls: f64, cad:String)  {
     // Create a drawing area
     let us1="BlandAltman".to_string() + &cad.to_string() +&"_".to_string() + &".png".to_string();

     let root = BitMapBackend::new(&us1, (800, 600)).into_drawing_area();
     root.fill(&WHITE).unwrap();
    
     let n=xx.len();
     let mut mixx=xx[0];
     let mut maxx=xx[0];
     let mut miyy=yy[0];
     let mut mayy=yy[0];
     for i in 0..n{
        if xx[i] < mixx {
            mixx=xx[i];
        }
        if xx[i] > maxx {
            maxx=xx[i];
        }
        if yy[i] < miyy {
            miyy =yy[i];
        }
        if yy[i] > mayy {
            mayy = yy[i];
        }
     }
     

     // Create a chart context
     let mut chart = ChartBuilder::on(&root)
         .caption("Bland - Altman", ("sans-serif", 40).into_font())
         .margin(10)
         .x_label_area_size(30)
         .y_label_area_size(30)
         .build_cartesian_2d(mixx..maxx, miyy..mayy)
         .unwrap();
 
     // Create a scatter series
     let mut data: Vec<(f64, f64)> = Vec::with_capacity (xx.len() as usize);
     for i in 0..xx.len() as usize{
        data.push((xx[i],yy[i]));
     }
     let mut data1: Vec<(f64, f64)> = Vec::with_capacity (xx.len() as usize);
     for i in 0..xx.len() as usize{
        data1.push((xx[i],med));
     }
     let mut data2: Vec<(f64, f64)> = Vec::with_capacity (xx.len() as usize);
     for i in 0..xx.len() as usize{
        data2.push((xx[i],ls));
     }
     let mut data3: Vec<(f64, f64)> = Vec::with_capacity (xx.len() as usize);
     for i in 0..xx.len() as usize{
        data3.push((xx[i],li));
     }
     // ejes
     let mut data4: Vec<(f64, f64)> = Vec::with_capacity (xx.len() as usize);
     for i in 0..xx.len() as usize{
        data4.push((xx[i],0.0));
     }
     let mut data5: Vec<(f64, f64)> = Vec::with_capacity (xx.len() as usize);
     for i in 0..xx.len() as usize{
        data5.push((0.0,yy[i]));
     }


     // Convert the data into a series of shape elements
     let scatter = data.into_iter().map(|(x, y)| Circle::new((x, y), 2, ORANGE.filled()));
 
     // Draw the scatter series on the chart
     chart.draw_series(scatter).unwrap();
    

     // Draw the average line
    chart.draw_series(std::iter::once(PathElement::new(
        data1,
        &BLUE,
    )))
    .unwrap();
     // Dibujar 
    chart.draw_series(std::iter::once(PathElement::new(
        data2,
        &RED,
    )))
    .unwrap();
    
    chart.draw_series(std::iter::once(PathElement::new(
        data3,
        &MAGENTA,
    )))
    .unwrap();

    chart.draw_series(std::iter::once(PathElement::new(
        data4,
        &BLACK,
    )))
    .unwrap();

    chart.draw_series(std::iter::once(PathElement::new(
        data5,
        &BLACK,
    )))
    .unwrap();

    }

fn normalk2(x:Vec<f64>) -> String{
        let resp1=" \n\n_Test de normalidad de D'Agostino - Pearson-k2 ";
        println!(" {} ",resp1.yellow());
        // el vector de datos tiene que estar ordenado
        let mut ord: Vec<f64>=x.clone();
        ord.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        let x=ord.clone();
        // inicia el análisis
        let alpha:f64 = 0.05;
        let n: f64 = ord.len() as f64;
        let s1: f64=suma(x.clone());
        let s2: f64=elevar_y_sumar(x.clone(), 2.0);
        let s3: f64 =elevar_y_sumar(x.clone(), 3.0);
        let s4: f64 =elevar_y_sumar(x.clone(), 4.0);
        let ss: f64 = s2-(s1.powf(2.0)/n);
        let v: f64 = ss/(n -1.0);
        let k3: f64 = ((n *s3)-(3.0*s1*s2)+((2.0*(s1.powf(3.0)))/n))/((n-1.0)*(n-2.0));
        //println!(" k3 {}",k3);
        let g1: f64 = k3/(v.powf(3.0)).powf(0.5);
        let k4: f64 = ((n +1.0)*((n * s4)-(4.0*s1*s3)+(6.0*(s1.powf(2.0))*(s2/n))-((3.0*(s1.powf(4.0)))/(n.powf(2.0))))/((n -1.0)*(n -2.0)*(n -3.0)))-((3.0*(ss.powf(2.0)))/((n -2.0)*(n -3.0)));
        let g2: f64 = k4/v.powf(2.0);
        let eg1: f64 = ((n-2.0)*g1)/(n*(n-1.0)).powf(0.5);
        //let eg2 = ((n-2.0)*(n-3.0)*g2)/((n+1.0)*(n-1.0))+((3.0*(n-1.0))/(n+1.0));
        let a: f64 = eg1* (((n+1.0)*(n+3.0))/(6.0*(n-2.0))).powf(0.5);
        let b: f64 = (3.0*((n.powf(2.0))+(27.0*n)-70.0)*((n+1.0)*(n+3.0)))/((n-2.0)*(n+5.0)*(n+7.0)*(n+9.0));
        let c: f64 = (2.0*(b-1.0)).powf(0.5)-1.0;
        let d: f64 = c.powf(0.5);
        let e: f64 = 1.0/((d.ln()).powf(0.5));
        let f: f64 = a/(2.0/(c-1.0)).powf(0.5);
        let zg1 = e*(f+(f.powf(2.0)+1.0).powf(0.5)).ln();
        println!(" zg1 {}",zg1);
        let g = (24.0*n*(n-2.0)*(n-3.0))/((n+1.0).powf(2.0)*(n+3.0)*(n+5.0));
        let h = ((n-2.0)*(n-3.0)*(g2).abs())/((n+1.0)*(n-1.0)*g.powf(0.5));
        let j = ((6.0*(n.powf(2.0)-(5.0*n)+2.0))/((n+7.0)*(n+9.0)))* ((6.0*(n+3.0)*(n+5.0))/((n*(n-2.0)*(n-3.0)))).powf(0.5);
        let k = 6.0+((8.0/j)*((2.0/j)+(1.0+(4.0/j.powf(2.0))).powf(0.5)));
        let l = (1.0 -(2.0/k))/(1.0 + h * (2.0/(k-4.0)).powf(0.5));
        let zg2 = (1.0 - (2.0/(9.0*k))-l.powf(1.0/3.0))/(2.0/(9.0*k)).powf(0.5);
        println!(" zg2 {}",zg2);
        let k2 = zg1.powf(2.0) + zg2.powf(2.0);
        let x2 = k2.clone();
        println!(" k2 {}",x2);

        let p: f64=1.0 - chi_ccdf(x2);
        println!(" P = {} ",p);
        if p>alpha{
            return "La muestra SI se distribuye normalmente ...".to_string()
            
        }else{
            return "La muestra NO viene de una distribución normal".to_string()
        }
        
}



fn mediana(data: &mut [f64]) -> f64 {
    
    if data.is_empty() {
        println!(" no hay Mediana ...");
    }
    
    data.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    
    let n = data.len();
    let mmed:f64;
    if n % 2 == 1 {
        mmed=data[n / 2];
    } else {
        mmed=(data[n / 2 - 1] + data[n / 2]) / 2.0;
    }
    mmed
}
fn quantile(v: &mut Vec<f64>, q: f64) -> f64 {
    // es correcta la petición?
    if q < 0.0 || q > 1.0 {
        println!(" No es correcto lo que pide ...");
    }
    // ordenar el vector
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    // encontrar indice y fracción
    let n = v.len() as f64;
    let i = (n - 1.0) * q;
    let f = i.fract();
    let k = i.floor() as usize;
    // interpolar el valor
    let fractil:f64;
    if k == n as usize - 1 {
        fractil=v[k];
    } else {
        fractil=v[k] + f * (v[k + 1] - v[k]);
    }
    fractil
}
fn pboxplot2(mat:&Vec<Vec<f64>>) {
    let resp1="_Detección de outliers ...";
    println!("\n\n{}",resp1.yellow());
    for columna in mat {
        // Convertir la columna en un vector
        let mut vector: Vec<f64> = columna.clone();

        // n
        let n=vector.len();
        // media
        let mut media=0.0;
        for i in 0..n{
            media=media+vector[i];
        }
        let _media=media/(n as f64);
        // mediana
        let _mediana=mediana(&mut vector);
        // cuartiles
        let q1=quantile(&mut vector, 0.25);
        let q2=quantile(&mut vector, 0.5);
        let q3=quantile(&mut vector, 0.75);
        let ri=q3-q1;
        // limites
        let ls=q3+1.5*ri;
        let li=q1-1.5*ri;
        // outliers
        let mut outliers: Vec<f64>=Vec::new();
        
        for i in 0..n{
            if vector[i] > ls || vector[i] < li{
                outliers.push(vector[i]);
            }
        }

        println!("  Outliers {:?}",outliers);
        println!("  Cuartiles {} {} {}  ls {} li {} \n\n",q1,q2,q3,ls,li);

        

    }
}

fn leeryasig(nom:String) -> (Vec<f64>, usize,f64,f64) {
    // Abre el archivo en modo de solo lectura (ignorando los errores).
    let file = File::open(nom.clone()).unwrap();
    let reader = BufReader::new(file);

    let mut v: Vec<f64> = Vec::new(); // crear vector para guardar la información

    println!("     ");
    // Lee el archivo línea por línea usando el iterador lines() de std::io::BufRead.
    for line in reader.lines() {
        // Convierte cada línea en un número flotante (ignorando los errores).
        let num: f64 = line.unwrap().parse().unwrap();
        // Guardar el número flotante.
        // println!("{}", num);
        v.push(num);

    }
    let resp1=" Datos ";
    let n=v.len();
    println!(" {} {:?} \n n= {}",resp1.yellow(),v,n);
    // guardar en matriz
   
    //minimo
    let mut minimo:f64=v[0];
    for i in 0..v.len() {
        if v[i] <= minimo {
            minimo=v[i];
        }
    }
    let resp1=" mínimo ";
    println!(" {} {:?}",resp1.yellow(),minimo);
    // maximo
    let mut maximo:f64=v[0];
    for i in 0..v.len() {
        if v[i] >= maximo {
            maximo=v[i];
        }
    }
    let resp1=" máximo ";
    println!(" {} {:?}",resp1.yellow(),maximo);
    /*
    Rango
    */
    let rango=maximo-minimo;
    let resp1=" rango ";
    println!(" {} {}",resp1.yellow(),rango);
    // aumentar
    let cinco=m_para_cinco(v.clone());
    // limites
    //println!(" cinco usado para ampliar limites {}",cinco);
    let cinco:f64 =5.0 * 10.0_f64.powf(-1.0 * (cinco as f64+1.0));
    let resp1=" cinco usado para ampliar limites ";
    println!(" {} {}",resp1.yellow(),cinco);
    let lsic=maximo + cinco as f64;
    let liic=minimo - cinco as f64;
    println!(" LIIC {} LSIC {}",liic,lsic);
    // numero de intervalos de clase
    
    let nn=f64::log10(n as f64);
    let ni = (3.322 * nn +1.0).round();
    let resp1=" número de intervalos de clase vía Sturges ";
    println!(" {} {}",resp1,ni);
    // amplitud de los intervalos de clase
    let aic=(lsic-liic)/ni;
    let resp1=" amplitud de los intervalos de clase ";
    println!(" {} {} ",resp1.yellow(),aic);
    let resp1="\n_Tabla de frecuencias ";
    println!(" {}",resp1.yellow());
    // intervalos de clase
    let mut l1 = Vec::new(); 
    let mut l2 = Vec::new(); 
    let mut pm=Vec::new(); 
    let mut frecc=Vec::new();
    let mut rel=Vec::new();
    let mut acu=Vec::new();


    l1.push(liic);
    for i in 1..(ni as usize +1){
        l2.push(liic + aic * i as f64);
        l1.push(liic + aic * i as f64);
    }
    let mut s:f64=0.0;
    let mut sacu:f64=0.0;
    for i in 0..ni as usize{
        pm.push((l1[i]+l2[i])/2.0); // puntos medios
        for j in 0..v.len() as usize{ // frecuencias
            if v[j]<=l2[i] && v[j]>l1[i]{
                s +=1.0;
            }    
        }
        frecc.push(s);
        rel.push(s/n as f64);
        sacu=sacu + s/n as f64;
        acu.push(sacu);
        s=0.0;
        //println!(" {} - {}  {}  {}",l1[i],l2[i],pm[i],frecc[i]);
    }


// Imprimir la matriz en forma de tabla
println!(" | {:^15} | {:^15} | {:^15} | {:^15} | {:^15} | {:^15}", "L1", "L2", "t","f","f*","F");
println!(" | {:^15} | {:^15} | {:^15} | {:^15} | {:^15} | {:^15}", "-", "-", "-","-","-","-");
for j in 0..ni as i32 {
    println!(" | {:^15.9} | {:^15.9} | {:^15.9} | {:^15} | {:^15.9} | {:^15.9}", l1[j as usize],l2[j as usize],pm[j as usize],frecc[j as usize],rel[j as usize],acu[j as usize]);
}
println!(" ");

let promedio:f64;
s=0.0;
for i in 0..=v.len() - 1{
    s=s+v[i];
}
promedio=s/v.len() as f64;
println!(" media  {}", promedio);
s=0.0;
for i in 0..=v.len() - 1{
    s=s+(v[i] - promedio).powf(2.0);
}
if n>=30 {
    s=s/n as f64;
}else if n<30 {
    s=s/(n as f64 -1.0);
}
let varianza:f64=s;
let ds:f64=s.powf(0.5);
println!(" var {}    ds {}",varianza,ds);

if n >= 8 {
    let nor =normalk2(v.clone());
    println!("  {}",nor);
}else {
    let resp1=" Se requiere muestra mayor e igual que 8 datos ... ";
    println!(" {}",resp1.red());
}

// tabla de valoración *************************************
// *********************************************************

let mut l1v:Vec<f64>=Vec::new();
let mut l2v:Vec<f64>=Vec::new();
let mut frec:Vec<f64>=Vec::new();

let ds: f64=varianza.powf(0.5);// desviación estandar

l1v.push(promedio+3.0*ds);
l1v.push(promedio+2.0*ds);
l1v.push(promedio+1.0*ds);
l1v.push(promedio);
l1v.push(promedio-1.0*ds);
l1v.push(promedio-2.0*ds);
l1v.push(promedio-3.0*ds);
l1v.push(-1e7);

l2v.push(1e7);
l2v.push(promedio+3.0*ds);
l2v.push(promedio+2.0*ds);
l2v.push(promedio+1.0*ds);
l2v.push(promedio);
l2v.push(promedio-1.0*ds);
l2v.push(promedio-2.0*ds);
l2v.push(promedio-3.0*ds);

let mut ss=0.0;
for j in 0..8{
    for i in 0..v.len(){
        if v[i]>l1v[j] && v[i]<=l2v[j] {
            ss +=1.0;
        }
    }
    frec.push(ss);
    ss=0.0;
}

let resp1="\n_Tabla de valoración";
println!(" {}",resp1.yellow());
println!(" | {:^15} | {:^15} | {:^15} | {:^15} |", "L1", "L2","frecuencia","Eval. Cualitativa");
for i in 0..8 {
    println!(" | {:^15.8} | {:^15.8} | {:^15.8} | ",l1v[i as usize],l2v[i as usize],frec[i as usize]);
}


(v,n,promedio,varianza)

}


fn cor(x:Vec<f64>, y:Vec<f64>)-> f64{
    let n=x.len();
    let mut s1=0.0;
    let mut s2=0.0;
    let mut s3=0.0;
    let mut s4=0.0;
    let mut s5=0.0;
    for i in 0..x.len(){
        s1=s1+ (x[i]*y[i]);         // s(x,y)
        s2=s2+ x[i].powf(2.0);      // s(x^2)
        s3=s3+ x[i];                // s(x)
        s4=s4+ y[i].powf(2.0);      // s(y^2)
        s5=s5+ y[i]                 // s(y)
    }
    let cor=(n as f64*s1-s3*s5)/((n as f64 * s2-s3.powf(2.0))*(n as f64 * s4-s5.powf(2.0))).powf(0.5);
    cor
}




use glob::glob;
use std::process::Command;


fn main(){
    let mut nombres: Vec<String> = Vec::new();  // lista de archivos de datos
    let mut matriz: Vec<Vec<f64>> = Vec::new();// matriz de esos datos
    let mut varis: Vec<f64> = Vec::new();// varianzas
    let mut medias: Vec<f64> = Vec::new();// promedios

    let mut enes:Vec<usize>=Vec::new();// tamaños de muestras

    // limpiar
    Command::new("cmd")
            .args(&["/C", "cls"])
            .status()
            .expect("Error al ejecutar el comando cls");
    // Iterar sobre los archivos que coinciden con el patrón
    for entrada in glob("numeros[0-9]*.txt").expect("Error al leer el patrón") {
        // Obtener el nombre del archivo como una cadena
        let nombre = entrada.unwrap().display().to_string();
        // Agregar el nombre al vector
        nombres.push(nombre);
    }

    // Mostrar el vector
    // println!("{:?}", nombres);
    let resp1="\n\n***** archivo ";
    for i in nombres.iter(){
        //leeryasig("numeros1.txt".to_string());
        println!("{} {} *****************************",resp1.blue(), i);
        let (vect,nn,prm,varianzas)=leeryasig(i.to_string());
        matriz.push(vect);                                    
        medias.push(prm);
        varis.push(varianzas);
        enes.push(nn);
    }

    // información para crear los boxplot de todas las muestras
    // estadisticos , no grafica
    pboxplot2(&matriz );
 
    let mut vee:Vec<f64>=Vec::new();
    for i in 0..enes[0]{
        vee.push(matriz[0][i]);
    }
    //println!("***** {:?}",vee);
    
    // verificar tamaños de muestras iguales
    let uu=enes[0];
    for i in 1..enes.len(){
        if uu != enes[i]{
            println!(" Tamaños de muestra diferentes ... ERROR ");
        }
    }
    let cols=nombres.len();
    let reng=uu;
    println!(" matriz :  [reng {} , cols {}]",reng, cols);

    // anova ******************************************************
    // ************************************************************

    // varianza total
    // media total
    let resp1: &str=" \n_Test ANOVA de una vía diseño completamente aleatorio";
    println!(" {}",resp1.yellow());
    varis.clear();
    let mut sum: f64=0.0;
    let mut sum2:f64=0.0;
    let mut n1=0.0;
    let mut n2=0.0;
    for i in 0..cols{
        for j in 0..reng{
            sum = sum + matriz[i][j];// suma de todos los datos de todas las muestras
            n1 +=1.0;
            sum2=sum2 + (matriz[i][j] - medias[i]).powf(2.0); // suma de cuadrados por muestra
            n2 +=1.0;
        }
        varis.push(sum2/(n2-1.));
        sum2=0.0;
        n2 =0.0;
    }
    let mediat: f64=sum/n1; // media total
    // varianza total
    sum=0.0;
    for i in 0..cols{
        for j in 0..reng{
            sum=sum+(matriz[i][j] - mediat).powf(2.0);// suma de cuadrados total
        }
    }
    let vart: f64=sum/(n1 -1.0);
    // mostrar medias y varianzas
    println!(" | {:^18} | {:^18} |", "Promedios","Varianzas");
    for i in 0..cols as usize{
        println!(" | {:^18.8} | {:^18.8} |", medias[i],varis[i]);
    }
    let sct=(n1 -1.0)*vart;
    sum=0.0;
    for i in 0..cols as usize{
        sum=sum + (enes[i] as f64 - 1.0)*varis[i];
    }
    let scd=sum;
    let sce=sct-scd;
    // fuente de variabilidad
    println!(" SCT {}\n SCD {}\n SCE {}",sct,scd,sce);
    // grados de libertad
    println!(" k-1= {} n-k= {} n-1= {}",cols-1,n1 as f64-cols as f64,n1 as f64-1.);
    // medios cuadrados
    let cme=sce/(cols as f64-1.0);
    let cmd=scd/(n1 as f64 - cols as f64);
    let f=cme/cmd;
    println!(" cme= {} cmd= {}  f= {}",cme,cmd,f);
    let df1=(cols as f64-1.0) as u64;
    let df2=(n1 as f64 - cols as f64) as u64;
    
    let p = 0.95;
    let q = f_quantile(p, df1, df2);
    println!(" EL {} percentil de la distribución F con g.l. ({}, {}) es {}", p * 100.0, df1, df2, q);
    if f > q {
        let resp1=" Al menos una media muestral es diferente...";
        println!(" {}",resp1.yellow());
        posthoc_scheffe(&medias,cmd,q);
    }else{
        let resp1=" Todas las medias muestrales son iguales...";
        println!(" {}",resp1.yellow());
        posthoc_scheffe(&medias,cmd,q);
    }
    
    // Test de Levene ***********************************************
    // **************************************************************
    // varianza total
    // media total
    let resp1=" \n_Test de Levene para verificar Homocedasticidad";
    println!(" {}",resp1.yellow());
    varis.clear();
    let matt: Vec<Vec<f64>>=matriz.clone(); // mantener datos originales
    for i in 0..cols{
        for j in 0..reng{
            matriz[i][j]=(matriz[i][j] - medias[i]).abs();
        }
    }

    let mut sum: f64=0.0;
    let mut sum2:f64=0.0;
    let mut n1=0.0;
    let mut n2=0.0;
    for i in 0..cols{
        for j in 0..reng{
            sum = sum + matriz[i][j];// suma de todos los datos de todas las muestras
            n1 +=1.0;
            sum2=sum2 + (matriz[i][j] - medias[i]).powf(2.0); // suma de cuadrados por muestra
            n2 +=1.0;
        }
        varis.push(sum2/(n2-1.));
        sum2=0.0;
        n2 =0.0;
    }
    let mediat: f64=sum/n1; // media total
    // varianza total
    sum=0.0;
    for i in 0..cols{
        for j in 0..reng{
            sum=sum+(matriz[i][j] - mediat).powf(2.0);// suma de cuadrados total
        }
    }
    let vart: f64=sum/(n1 -1.0);
    // mostrar medias y varianzas
    println!(" | {:^18} | {:^18} |", "Promedios","Varianzas");
    for i in 0..cols as usize{
        println!(" | {:^18.8} | {:^18.8} |", medias[i],varis[i]);
    }
    let sct=(n1 -1.0)*vart;
    sum=0.0;
    for i in 0..cols as usize{
        sum=sum + (enes[i] as f64 - 1.0)*varis[i];
    }
    let scd=sum;
    let sce=sct-scd;
    // fuente de variabilidad
    println!(" SCT {}\n SCD {}\n SCE {}",sct,scd,sce);
    // grados de libertad
    println!(" k-1= {} n-k= {} n-1= {}",cols-1,n1 as f64-cols as f64,n1 as f64-1.);
    // medios cuadrados
    let cme=sce/(cols as f64-1.0);
    let cmd=scd/(n1 as f64 - cols as f64);
    let f=cme/cmd;
    println!(" cme= {} cmd= {}  f= {}",cme,cmd,f);
    let df1=(cols as f64-1.0) as u64;
    let df2=(n1 as f64 - cols as f64) as u64;
    
    let p = 0.95;
    let q = f_quantile(p, df1, df2);
    println!(" EL {} percentil de la distribución F con g.l. ({}, {}) es {}", p * 100.0, df1, df2, q);
    if f > q {
        let resp1: &str=" Al menos una varianza muestral es diferente...";
        println!(" {}",resp1.yellow());
    }else{
        let resp1: &str=" Todas las varianzas muestrales son iguales...";
        println!(" {}",resp1.yellow());
    }

    // crear todas las combinaciones de comparaciones por parejas
    // para los casos de correlación y coeficiente de concordancia de
    // Lawrence - Linn
    matriz.clear();
    let mut numerosc:Vec<usize>=Vec::new();
    let mut vec11:Vec<f64>=Vec::new();
    let mut vec22:Vec<f64>=Vec::new();
    for i in 0..cols{
        numerosc.push(i);
    }
    let resp1=" \n_Combinaciones de las muestras";
    println!(" {} {:?}",resp1.yellow(),numerosc);
    for i in 0..cols - 1 {
        let a = numerosc[i];
        for j in i + 1..cols {
            let b = numerosc[j];
            // combinación de a y b
            println!(" {} vs {}", a, b);
            let cad=a.to_string()+ &"_".to_string()+&b.to_string();
            for z in 0..reng{
                vec11.push(matt[a][z]);
                vec22.push(matt[b][z]);
            }
            let c=cor(vec11.clone(),vec22.clone());
            let resp1=" Correlación de Pearson ";
            println!(" {} {}",resp1.yellow(),c);
            let ccc = l_linn(&vec11, &vec22);
            let resp1=" Concordancia de Lawrence - Linn ";
            println!(" {}  {} ",resp1.yellow(),ccc);
            // datos de Blant - Altman
            blanalt(&vec11,&vec22,cad);
            vec11.clear();
            vec22.clear();
        }
    }

    // Mantener la consola abierta después de salir
    Command::new("cmd")
            .args(&["/C", "cmd /k"])
            .status()
            .expect("Error al ejecutar el comando cmd /k");

    
}
