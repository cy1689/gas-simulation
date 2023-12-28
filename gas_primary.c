/**********************************************************************/
/*                                                                    */
/* Compile this file using 'gcc gas_primary.c -lm -o gp'.             */
/* Then run './gp' to generate all the data presented in the papers.  */
/*                                                                    */
/* When the simulation time is set to 6.0e8, it may take a long time. */
/* You can modify this to a smaller value, for example, 6e5. The      */
/* results won't be as accurate, but it will produce all the results  */
/* in a shorter time.                                                 */
/*                                                                    */
/**********************************************************************/


#define SIMULATION_TIME  6.0e8

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <float.h>

#ifndef max
#define max(A, B) ((A)>(B) ? (A) :(B))
#endif

#ifndef min
#define min(A, B) ((A)<(B) ? (A) :(B))
#endif

#define TOP (-1)
#define RIGHT (-2)
#define BOTTOM (-3)
#define LEFT (-4)
#define FRONT (-5)
#define BACK (-6)

#define SQRT3 1.7320508075689

typedef struct {
  /* for testing */
  int enable_particle_particle_collision;
  int is_3D;
  int do_density_stat;
  int density_at_once;
  int is_interior_p;
  double radius;
  double mass;
  /* this is with respect to the reduced container corner (2D) or edge (3D) */
  double excluded_d;
  double density_prediction;
} Setting;

typedef struct {
  double startTime;
  double startX;
  double startY;
  double startZ; 
  double velocityX;
  double velocityY;
  double velocityZ;
  double collisionTime;
  int collideWith;
} Particle;

typedef struct {
  double x1, x2, y1, y2, z1, z2;
  double volumeRatio;
  double inTime;
  double inTimeVySq; 
} PositionTime;

typedef struct {
  double w, h, l;
} Container;

typedef struct {
  double x, y, z;
} Coordinate;

typedef struct {
  double impulseX, impulseY, impulseZ;
} Impulse;

typedef struct {
  double c0, c1, c2;
  double rmse;
  double c0Sd;
} QuadraticFittingCoef;

double myRand() {
  return (double) random()/(1L<<31);
}

void averageAndSd(double *ave, double *var, double *buf, int n) {
  double sum=0., sum2=0.;
  int i;

  for(i=0; i<n; i++) {
    sum += buf[i];
    sum2 += buf[i]*buf[i];
  }

  sum /= (double) n;
  sum2 /= (double) n;

  *ave = sum;
  // roundoff error may result in sum2 being less than sum*sum
  sum2 -= sum*sum;
  *var = sum2 <=0. ? 0. :sqrt(sum2/(double) (n-1));
}

double matrix3x3Det (
  double a11,
  double a12,
  double a13,
  double a21,
  double a22,
  double a23,
  double a31,
  double a32,
  double a33) {

  return a11*a22*a33 + a12*a23*a31 + a13*a32*a21
    - a13*a22*a31 - a23*a32*a11 - a33*a21*a12;
}

void quadraticFitting(QuadraticFittingCoef *coef,
  double *x,
  double *y,
  double *yVar,
  int n) {
  int i;
  double xSum, x2Sum, x3Sum, x4Sum;
  double ySum, yxSum, yx2Sum;
  double A, B, C, D;

  xSum = 0.;
  x2Sum = 0.;
  x3Sum = 0.;
  x4Sum = 0.;
  ySum = 0.;
  yxSum = 0.;
  yx2Sum = 0.;

  for(i=0; i<n; i++) {
    xSum += x[i];
    x2Sum += x[i]*x[i];
    x3Sum += x[i]*x[i]*x[i];
    x4Sum += x[i]*x[i]*x[i]*x[i];
    ySum += y[i];
    yxSum += y[i]*x[i];
    yx2Sum += y[i]*x[i]*x[i];
  }

  D = matrix3x3Det ((double) n, xSum, x2Sum,
    xSum, x2Sum, x3Sum,
    x2Sum, x3Sum, x4Sum);
  A = matrix3x3Det (ySum, xSum, x2Sum,
    yxSum, x2Sum, x3Sum,
    yx2Sum, x3Sum, x4Sum);
  B = matrix3x3Det ((double) n, ySum, x2Sum,
    xSum, yxSum, x3Sum,
    x2Sum, yx2Sum, x4Sum);
  C = matrix3x3Det ((double) n, xSum, ySum,
    xSum, x2Sum, yxSum,
    x2Sum, x3Sum, yx2Sum);

  coef->c0 = A/D;
  coef->c1 = B/D;
  coef->c2 = C/D;

  coef->rmse = 0.;
  for(i=0; i<n; i++) {
    double t = coef->c0+coef->c1*x[i] + coef->c2*x[i]*x[i]-y[i];
    coef->rmse += t*t;
  }
  coef->rmse /= n-3;
  coef->rmse = sqrt(coef->rmse);

  coef->c0Sd = 0.;
  if(yVar) {
    for(i=0; i<n; i++) {
      double E = matrix3x3Det (1, xSum, x2Sum,
        x[i], x2Sum, x3Sum,
        x[i]*x[i], x3Sum, x4Sum);
      coef->c0Sd += E*E*yVar[i];
    }
    coef->c0Sd = sqrt(coef->c0Sd) / D;
  }
}

double particleCollisionTime(Setting *setting, Particle *a, Particle *b) {
  double startT;
  double dta, xa, ya, za;
  double dtb, xb, yb, zb;
  double dx, dy, dz;
  double dvx, dvy, dvz;
  double dv;

  startT = a->startTime;
  if(b->startTime > a->startTime) {
    startT = b->startTime;
  }

  dta = startT - a->startTime;
  xa = a->startX + dta * a->velocityX;
  ya = a->startY + dta * a->velocityY;
  if(setting->is_3D) {
    za = a->startZ + dta * a->velocityZ;
  }

  dtb = startT - b->startTime;
  xb = b->startX + dtb * b->velocityX;
  yb = b->startY + dtb * b->velocityY;
  if(setting->is_3D) {
    zb = b->startZ + dtb * b->velocityZ;
  }

  // b's displacement with respect to a
  dx = xb - xa;
  dy = yb - ya;
  if(setting->is_3D) {
    dz = zb - za;
  }

  // b's velocity with respect to a
  dvx = b->velocityX - a->velocityX;
  dvy = b->velocityY - a->velocityY;
  if(setting->is_3D) {
    dvz = b->velocityZ - a->velocityZ;

    dv = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
  } else {
    dv = sqrt(dvx*dvx + dvy*dvy);
  }
 
  if(dv == 0.) {
    return DBL_MAX;
  }

  // two particles are moving toward each other only if the dot 
  // product of their velocity and displacement is negative.
  double dotP = dx * dvx + dy * dvy;
  if(setting->is_3D) {
    dotP += dz * dvz;
  }

  if(dotP >= 0.) {
    return DBL_MAX;
  }

  double dis;
  if(setting->is_3D) {
    dis = sqrt(dx*dx + dy*dy + dz*dz);
  } else {
    dis = sqrt(dx*dx + dy*dy);
  }

  // distance in the direction of velocity 
  double vDis = -dotP/dv;
  double t = dis*dis -vDis*vDis;
  double vPerpendicularDis = t<=0. ? 0. :sqrt(t);
  double tt = 4*setting->radius*setting->radius -
    vPerpendicularDis*vPerpendicularDis;
  if(tt<=0.) {
    return DBL_MAX;
  }

  return startT + (vDis - sqrt(tt))/dv;
}

void calcCollisionTime(Setting *setting,
  Particle *p, int n, int idx, Container *container)
{
  double t = DBL_MAX;
  double tt;
  int collideWith = TOP;
  int i;

  if(p[idx].velocityX > 0) {
    tt = (container->w- p[idx].startX - setting->radius)/p[idx].velocityX;
    if(tt<t) {
      t = tt;
      collideWith = RIGHT;
    }
  }
  else if(p[idx].velocityX < 0) {
    tt = (setting->radius - p[idx].startX)/p[idx].velocityX;
    if(tt<t) {
      t = tt;
      collideWith = LEFT;
    }
  }

  if(p[idx].velocityY > 0) {
    tt = (container->h- p[idx].startY - setting->radius)/p[idx].velocityY;
    if(tt<t) {
      t = tt;
      collideWith = BOTTOM;
    }
  }
  else if(p[idx].velocityY < 0) {
    tt = (setting->radius - p[idx].startY)/p[idx].velocityY;
    if(tt<t) {
      t = tt;
      collideWith = TOP;
    }
  }

  if(setting->is_3D) {
    if(p[idx].velocityZ > 0) {
      tt = (container->l- p[idx].startZ - setting->radius)/p[idx].velocityZ;
      if(tt<t) {
        t = tt;
        collideWith = BACK;
      }
    }
    else if(p[idx].velocityZ < 0) {
      tt = (setting->radius - p[idx].startZ)/p[idx].velocityZ;
      if(tt<t) {
        t = tt;
        collideWith = FRONT;
      }
    }
  }

  if(t != DBL_MAX) {
    t += p[idx].startTime;
  }

  if(setting->enable_particle_particle_collision){
    for(i=0; i<n; i++) {
      if(i != idx) {
        tt = particleCollisionTime(setting, p+idx, p+i);
        if(tt < t) {
          t = tt;
          collideWith = i;
        }
      }
    }
  }

  p[idx].collideWith = collideWith;
  p[idx].collisionTime = t;
}

int overlap(Setting *setting, Particle *a, Particle *b) {
  double dx = b->startX - a->startX;
  double dy = b->startY - a->startY;
  double sum = dx*dx + dy*dy;
  double dz;

  if(setting->is_3D) {
    dz = b->startZ - a->startZ;
    sum += dz*dz;
  }

  return (sum < setting->radius*setting->radius *4.0 ? 1:0);
}

int overlapAny(Setting *setting, Particle *arr, int n, Particle *a) {
  int i;

  for(i=0; i<n; i++) {
    if(overlap(setting, arr+i, a)) {
      return 1;
    }
  }

  return 0;
}

double volume_0(Setting *setting, Container *container) {
  double r = setting->radius;
  double v = (container->w - r*2.0)*(container->h - r*2.0);

  if(setting->is_3D) {
    v *= container->l - r*2.0;
  }

  return v;
}

double pressure_0(Setting *setting, Container *container, int particleNum) {
  double v0 = volume_0(setting, container);
  double p0;

  // average v^2 is 1
  if(setting->is_3D) {
    p0 = particleNum * setting->mass * (1.0/3.0)/v0;
  } else {
    p0 = particleNum * setting->mass * 0.5/v0;
  }
  return p0;
}

double impulsePerUnitArea(Setting *setting, Impulse *impulse, Container *container) {
  double r = setting->radius;
  double totImpulse, totPressure, totSurf;
  double fac;

  fac = setting->is_interior_p ? 2.0*setting->excluded_d : 0.0;
  fac += 2.0;

  if(setting->is_3D) {
    totPressure  = impulse->impulseX /(2.0*(container->h - r*fac)*(container->l - r*fac));
    totPressure += impulse->impulseY /(2.0*(container->w - r*fac)*(container->l - r*fac));
    totPressure += impulse->impulseZ /(2.0*(container->w - r*fac)*(container->h - r*fac));
    return totPressure/3.0;
  } 

  totPressure  = impulse->impulseX /(2.0*(container->h - r*fac));
  totPressure += impulse->impulseY /(2.0*(container->w - r*fac));
  return totPressure * 0.5;
}

Particle *allocateParticles(int n) {
  return (Particle *) malloc(sizeof(Particle)*n);
}

void initParticles(Setting *setting, Particle *p, int n, Container *container)
{
  int i;
  double angle;
  double z;

  for (i=0; i<n; i++) {
    p[i].startTime = 0.;
    angle = M_PI*2.0 * myRand();
    p[i].velocityX = cos(angle);
    p[i].velocityY = sin(angle);
    //printf("vy^2 = %f\n", p[i].velocityY*p[i].velocityY);
    if(setting->is_3D) {
      z = 2.0 * (myRand()-0.5f);
      p[i].velocityZ = z;
      double t = 1.0-z*z;
      // t should be >= 0, but roundoff error may make it <0
      z = t<=0. ? 0. : sqrt(t);
      p[i].velocityX *= z;
      p[i].velocityY *= z;
    }

    do {
      p[i].startX = setting->radius + (container->w-2*setting->radius) * myRand();
      p[i].startY = setting->radius + (container->h-2*setting->radius) * myRand();
      if(setting->is_3D) {
        p[i].startZ = setting->radius + (container->l-2*setting->radius) * myRand();
      }
    } while(overlapAny(setting, p, i, p+i));
  }
}

int nextCollisionParticle(Particle *p, int n) {
  double collisionTime = p[0].collisionTime;
  int i, idx = 0;

  for(i=1; i<n; i++) {
    if(p[i].collisionTime < collisionTime) {
      idx = i;
      collisionTime = p[i].collisionTime;
    }
  }

  return idx;
}

PositionTime *allocateDensityStatisticsBoxes(int n) {
  PositionTime *pt;
  pt = (PositionTime *) malloc(sizeof(PositionTime)*n);
  if(!pt) {
    fprintf(stderr, "Failed to allocate density statistics boxes\n");
  }

  return pt;
}

void printDensityStatisticsBoxesInfo(Setting *setting, PositionTime *pt, int cnt, Container *container) {
  printf("Container w=%f h=%f", container->w, container->h);
  if(setting->is_3D) {
    printf(" l=%f", container->l);
  }

  printf(". Particle radius: %f\n", setting->radius);
  printf("Reduced container w=%f h=%f", 
    container->w-2*setting->radius,
    container->h-2*setting->radius);
  if(setting->is_3D) {
    printf(" l=%f", container->l-2*setting->radius);
  }
  printf("\n");

  for (int i=0; i<cnt; i++) {
    printf("Box %d info%s:\n",
      i,
      i==0 ? " (coordinates are with respect to the reduced container)":""
    );
    printf("    x1=%f x2=%f\n", 
       pt[i].x1-setting->radius, pt[i].x2-setting->radius);
    printf("    y1=%f y2=%f\n", 
       pt[i].y1-setting->radius, pt[i].y2-setting->radius);
    if(setting->is_3D) {
      printf("    z1=%f z2=%f\n", 
        pt[i].z1-setting->radius, pt[i].z2-setting->radius);
    }
    //printf("    %s with respect to reduced container %f\n", 
    //  setting->is_3D ? "Volume" : "Area",
    //  pt[i].volumeRatio);
  }
}


void initDensityStatisticsBoxes(Setting *setting, 
  PositionTime *pt, 
  int cnt, 
  double *ptThicknessBuf,
  Container *container, 
  int particleNum) {
  double r = setting->radius;
  double w = container->w, h=container->h, l, limit; 
  double v0 = volume_0(setting, container);
  int i;

  if(setting->is_3D) {
    l = container->l;
  }

  limit = r + setting->excluded_d;

  for(i=0; i<cnt; i++) {
    pt[i].x1 = limit;
    pt[i].x2 = w-limit;
    pt[i].y1 = setting->radius; 
    pt[i].y2 = pt[i].y1 + setting->radius*0.04 * (1.0-(double) i/cnt);
    pt[i].inTime = 0.;
    pt[i].inTimeVySq = 0.;
    pt[i].volumeRatio = (pt[i].x2-pt[i].x1)*(pt[i].y2-pt[i].y1)/v0;
    if(setting->is_3D) {
      pt[i].z1 = limit;
      pt[i].z2 = l-limit;
      pt[i].volumeRatio *= pt[i].z2-pt[i].z1;
    }
    ptThicknessBuf[i] = pt[i].y2 - pt[i].y1;
  }
  printDensityStatisticsBoxesInfo(setting, pt, cnt, container);
}

void initSetting(Setting *setting) {
  setting->enable_particle_particle_collision = 1;
  setting->is_3D = 0;
  setting->do_density_stat = 0;
  setting->density_at_once = 0;
  setting->is_interior_p = 0;
  setting->radius = 1.0;
  setting->mass = 1.0;
  setting->excluded_d = 0.0;
}

void initImpulse(Impulse *impulse) {
  impulse->impulseX = 0.;
  impulse->impulseY = 0.;
  impulse->impulseZ = 0.;
}

void freeInTimeBuf(double **buf, int positionTimeCount) {
  int i;

  if(!buf) return;

  for(i=0; i<positionTimeCount; i++) {
    free(buf[i]);
  }
  free(buf);
}

double **allocateInTimeBuf(int positionTimeCount, int particleSimulationCount) {
  double **inTimeBuf;
  int i;

  inTimeBuf = (double **) malloc(sizeof(double *)*2*positionTimeCount);
  if(!inTimeBuf) {
    fprintf(stderr, "Failed to allocate inTimeBuf\n");
    return 0;
  }

  for(i=0; i < 2*positionTimeCount; i++) {
    inTimeBuf[i] = (double *) malloc(sizeof(double)*particleSimulationCount);
    if(!inTimeBuf[i]) {
      fprintf(stderr, "Failed to allocate inTimeBuf[]\n");
      freeInTimeBuf(inTimeBuf, i);
      return 0;
    }
  }

  return inTimeBuf;
}

int intersectT(double *t1x,
  double *t2x,
  double ptx1,
  double ptx2,
  double x1,
  double x2) {
  double ptxFirst, ptxSecond;

  if(x1==x2) {
    if(x1<ptx1 || x1 >= ptx2) {
      return 0;
    }

    *t1x = 0.;
    *t2x = 1.;

    return 1;
  }

  if(x1 < x2) {
    ptxFirst = ptx1;
    ptxSecond = ptx2;
  } else {
    ptxFirst = ptx2;
    ptxSecond = ptx1;
  }

  *t1x = (ptxFirst-x1)/(x2-x1);
  *t2x = (ptxSecond-x1)/(x2-x1);

  return *t1x<1. && *t2x > 0.;
}

void updatePositionTime(Setting *setting,
  PositionTime *pt,
  int n,
  Particle *p,
  double dt) {
  double x1 = p->startX;
  double y1 = p->startY;
  double x2 = x1 + dt * p->velocityX;
  double y2 = y1 + dt * p->velocityY;
  double z1, z2;
  double t1x, t2x;
  double t1y, t2y;
  double t1z, t2z;
  double t1, t2;
  int i;

  if(setting->is_3D) {
    z1 = p->startZ;
    z2 = z1 + dt * p->velocityZ;
  }

  for(i=0; i<n; i++) {
    // without this, the code still works. When the sampling region is
    // very thin in y direction. Most trajectories satisfy the following,
    // so this speeds up the code.
    if(y1<=y2 ? y2 <= pt[i].y1 : y2 >= pt[i].y2) {
      continue;
    }

    if(!intersectT(&t1x, &t2x, pt[i].x1, pt[i].x2, x1, x2)) {
      continue;
    }

    t1 = max(0., t1x);
    t2 = min(1., t2x);

    if(!intersectT(&t1y, &t2y, pt[i].y1, pt[i].y2, y1, y2)) {
      continue;
    }

    t1 = max(t1, t1y);
    t2 = min(t2, t2y);
    if(t1 >= t2) {
      continue;
    }

    if(setting->is_3D) {
      if(!intersectT(&t1z, &t2z, pt[i].z1, pt[i].z2, z1, z2)) {
        continue;
      }

      t1 = max(t1, t1z);
      t2 = min(t2, t2z);
      if(t1 >= t2) {
        continue;
      }
    }

    pt[i].inTime += dt * (t2-t1);
    pt[i].inTimeVySq += dt * (t2-t1) * p->velocityY * p->velocityY;
  } 
}

void updatePositionTimeLast(Setting *setting,
  PositionTime *pt,
  int positionTimeCount,
  Particle *p,
  int particleNum,
  double simulationT) {
  double dt;
  int i;

  for(i=0; i<particleNum; i++) {
    dt = simulationT - p[i].startTime;

    if(dt > 0.) {
      updatePositionTime(setting, pt, positionTimeCount, p + i, dt);
    }
  }
}

void processCollision(Impulse *impulse,
  Setting *setting,
  Particle *p, int n, int idx, Container *container,
  PositionTime *pt, int positionTimeCount
) {
  double collisionTime = p[idx].collisionTime;
  int i = p[idx].collideWith;
  int j;
  double dt;
  double r = setting->radius; 
  double excluded_d = setting->excluded_d + r;
  double dx, dy, dz, ds2;
  double dvx, dvy, dvz, dotP;
  double factor;

  dt = collisionTime - p[idx].startTime;
  if(setting->do_density_stat) {
    updatePositionTime(setting, pt, positionTimeCount, p+idx, dt);
  }

  p[idx].startTime = collisionTime;
  p[idx].startX += dt * p[idx].velocityX;
  p[idx].startY += dt * p[idx].velocityY;
  if(setting->is_3D) {
    p[idx].startZ += dt * p[idx].velocityZ;
  }

  if(i<0) {
    if(i==LEFT || i==RIGHT) {
      if(impulse && (!setting->is_interior_p ||
         (p[idx].startY >= excluded_d && p[idx].startY < container->h-excluded_d
        && (!setting->is_3D ||
          (p[idx].startZ >= excluded_d && p[idx].startZ < container->l-excluded_d))))){ 
        impulse->impulseX += fabs(2.0*setting->mass*p[idx].velocityX);
      }
      p[idx].velocityX *= -1.0;
    } else if(i==TOP || i==BOTTOM) {
      if(impulse && (!setting->is_interior_p ||
        (p[idx].startX >= excluded_d && p[idx].startX < container->w-excluded_d
        && (!setting->is_3D || 
          (p[idx].startZ >= excluded_d && p[idx].startZ < container->l-excluded_d))))){ 
        impulse->impulseY += fabs(2.0*setting->mass*p[idx].velocityY);
      }

      p[idx].velocityY *= -1.0;
    } else {
      if(i != FRONT && i != BACK) {
        fprintf(stderr, "Collision with %d\n", i);
        exit(1);
      }

      if(impulse && (!setting->is_interior_p ||
        (p[idx].startX >= excluded_d && p[idx].startX < container->w-excluded_d && 
        p[idx].startY >= excluded_d && p[idx].startY < container->h-excluded_d))) { 
        impulse->impulseZ += fabs(2.0*setting->mass*p[idx].velocityZ);
      }

      p[idx].velocityZ *= -1.0;
    } 
  } else {
    dt = collisionTime - p[i].startTime;
    if(setting->do_density_stat) {
      updatePositionTime(setting, pt, positionTimeCount, p+i, dt);
    }
    p[i].startTime = collisionTime;
    p[i].startX += dt * p[i].velocityX;
    p[i].startY += dt * p[i].velocityY;
    if(setting->is_3D) {
      p[i].startZ += dt * p[i].velocityZ;
    }

    // use formula from https://en.wikipedia.org/wiki/Elastic_collision
    dx = p[i].startX-p[idx].startX;
    dy = p[i].startY-p[idx].startY;
    if(setting->is_3D) {
      dz = p[i].startZ-p[idx].startZ;
    }
    ds2 = dx*dx + dy*dy;
    if(setting->is_3D) {
      ds2 += dz*dz;
    }

    dvx = p[i].velocityX - p[idx].velocityX;
    dvy = p[i].velocityY - p[idx].velocityY;
    if(setting->is_3D) {
      dvz = p[i].velocityZ - p[idx].velocityZ;
    }

    dotP = dx*dvx + dy*dvy;
    if(setting->is_3D) {
      dotP += dz*dvz;
    }

    factor = 2*dotP / (setting->mass * 2 * ds2);
    p[idx].velocityX += factor * setting->mass * dx;
    p[idx].velocityY += factor * setting->mass * dy;
    if(setting->is_3D) {
      p[idx].velocityZ += factor * setting->mass * dz;
    }

    p[i].velocityX -= factor * setting->mass * dx;
    p[i].velocityY -= factor * setting->mass * dy;
    if(setting->is_3D) {
      p[i].velocityZ -= factor * setting->mass * dz;
    }
  }

  for(j=0; j<n; j++) {
    // i can be a wall
    if(j==idx || j==i ||
      p[j].collideWith == idx || (i>=0 && p[j].collideWith==i))
    {
      calcCollisionTime(setting, p, n, j, container);
    }
  }
}

void oneSimulation(Impulse *totImpulse,
  Setting *setting, 
  Particle *p, int particleNum, Container *container,
  double simulationT,
  PositionTime *pt, int positionTimeCount
) {
  double t;
  int i;

  if(setting->do_density_stat) {
    for(i=0; i<positionTimeCount; i++) {
      pt[i].inTime = 0.;
      pt[i].inTimeVySq = 0.;
    }
  }

  initParticles(setting, p, particleNum, container);
  for(i=0; i<particleNum; i++) {
    calcCollisionTime(setting, p, particleNum, i, container);
  }

  while(1) {
    i = nextCollisionParticle(p, particleNum);
    t = p[i].collisionTime;
    if(t > simulationT) {
      if(setting->do_density_stat) {
        updatePositionTimeLast(setting, pt, positionTimeCount, p, 
        particleNum, simulationT);
      }
      break;
    }

    processCollision(totImpulse, setting, p, particleNum, i, container,
      pt, positionTimeCount
    );
  }
}

double alpha2D_12(double w_2r, double h_2r, double r) {
  // note it's r (instead of 2r) in the following power terms
  double r_2 = r*r;
  double r_3 = r_2*r;
  double r_4 = r_3*r;

  double v0 = w_2r*h_2r;

  return (4.0*M_PI * r_2 - 32.0*((w_2r+h_2r)/3.0*r_3 - r_4/4.0)/v0)/v0;
}

double beta2D_12_at_x1(Setting *setting, double w_2r, double h_2r, double r) {
  if(setting->is_interior_p) {
    return 2.0*M_PI/(w_2r*h_2r)*r*r;
  }

  // note it's r (instead of 2r) in the following
  return (2.0*M_PI * h_2r - 16.0/3.0*r)*r*r/w_2r/h_2r/h_2r; 
}

double alpha2D(int pN, double w_2r, double h_2r, double r) {
  // for pN particles, there are 2 out of pN pairs, that is 
  // cN = pN*(pN-1)/2. Each pair has one primary inaccessible region.
  double tot = (double) (pN*(pN-1)/2) * alpha2D_12(w_2r, h_2r, r);
  //printf("alpha 0: %f\n", tot);

  return tot;
}

// note container surface "area" for 2D is actually surface length.
// pN is particle number
double beta2D_at_x(Setting *setting, int pN, double w_2r, double h_2r, double r) {
  if(!setting->is_interior_p && pN > 2) {
    fprintf(stderr, "When the number of particles exceeds two, the pressure is calculated only for the interior boundary region.\n");
    exit(1);
  }

  double tot = alpha2D(pN, w_2r, h_2r, r);
  //printf("0x: %f\n", tot);
  tot += (double) (pN-1) *
    (beta2D_12_at_x1(setting, w_2r, h_2r, r)-alpha2D_12(w_2r, h_2r, r));

  //printf("1x: %f\n", tot);
  
  return tot;
}

double alpha3D_12(double w_2r, double h_2r, double l_2r, double r) {
  double v0 = w_2r*h_2r*l_2r;

  double r_3 = r*r*r;
  double r_4 = r_3*r;
  double r_5 = r_4*r;
  double r_6 = r_5*r;
  return (32.0*M_PI/3.0*r_3-M_PI*8.0*(1.0/w_2r+1.0/h_2r+1.0/l_2r)*r_4 +
    (256.0/15.0*(w_2r+h_2r+l_2r)*r_5-32.0/3*r_6)/v0)/v0;
}

double beta3D_12_at_x1(Setting *setting, double w_2r, double h_2r, double l_2r, double r) {
  double v0 = w_2r*h_2r*l_2r;

  if(setting->is_interior_p) {
    return (M_PI*16.0/3.0/v0)*r*r*r;
  }

  return (M_PI*16.0/3.0-4.0*(M_PI*(h_2r+l_2r)*r-32.0/15.0*r*r)/(h_2r*l_2r))/v0*r*r*r;
}

double alpha3D(int pN, double w_2r, double h_2r, double l_2r, double r) {
  // for pN particles, there are 2 out of pN pairs, that is 
  // cN = pN*(pN-1)/2. Each pair corresponds to one primary region.
  double tot = (double) (pN*(pN-1)/2) * alpha3D_12(w_2r, h_2r, l_2r, r);

  return tot;
}

double beta3D_at_x(Setting *setting, int pN, double w_2r, double h_2r, double l_2r, double r) {
  if(!setting->is_interior_p && pN > 2) {
    fprintf(stderr, "For particle number greater than 2, "
      "is_interior_p must be true\n");
    exit(1);
  }

  double tot = alpha3D(pN, w_2r, h_2r, l_2r, r);
  tot += (double) (pN-1) *
    (beta3D_12_at_x1(setting, w_2r, h_2r, l_2r, r)-alpha3D_12(w_2r, h_2r, l_2r, r));

  return tot;
}

double analyticalBoundaryDensityPrimary(Setting *setting, 
  int particleNum, 
  Container *container) {
  double r = setting->radius;
  double w_2r = container->w-2.0*r; 
  double h_2r = container->h-2.0*r;
  double ratio = (2.0*r)*(2.0*r)/(w_2r*h_2r);
  double alpha_12, alpha;

  if(setting->is_3D) {
    double l_2r = container->l - 2.0*r;
    alpha_12 = alpha3D_12(w_2r, h_2r, l_2r, r);
    alpha = (double) particleNum*(particleNum-1)/2 * alpha_12;
    //printf("alpha_12=%f alpha=%f\n", alpha_12, alpha);

    ratio *= 2.0*r/l_2r;

    return (1.0-alpha-(M_PI*2.0/3.0*ratio-alpha_12)*(particleNum-1))/
      (1.0-alpha);
  }

  alpha_12 = alpha2D_12(w_2r, h_2r, r);
  alpha = (double) particleNum*(particleNum-1)/2 * alpha_12;

  return (1.0-alpha-(M_PI/2.0*ratio-alpha_12)*(particleNum-1))/
         (1.0-alpha);
}

double analyticalBoundaryDensityPrimaryS(Setting *setting, 
  int particleNum, 
  Container *container) {
  double r = setting->radius;
  double w_2r = container->w-2.0*r; 
  double h_2r = container->h-2.0*r;
  double ratio = (2.0*r)*(2.0*r)/(w_2r*h_2r);
  double alpha_12;

  if(setting->is_3D) {
    double l_2r = container->l - 2.0*r;

    ratio *= 2.0*r/l_2r;

    alpha_12 = alpha3D_12(w_2r, h_2r, l_2r, r);
    return 1.0-(M_PI*2.0/3.0*ratio-alpha_12)*(particleNum-1);
  }

  alpha_12 = alpha2D_12(w_2r, h_2r, r);
  return 1.0-(M_PI/2.0*ratio-alpha_12)*(particleNum-1);
}

double analyticalPressure(Setting *setting, int particleNum, Container *container) {
  double analyticalP;

  if(setting->is_interior_p) {
    analyticalP = particleNum > 2
        ? analyticalBoundaryDensityPrimaryS(setting, particleNum, container)
        : analyticalBoundaryDensityPrimary(setting, particleNum, container);
  } else {
    double r = setting->radius;
    double w_2r = container->w - r*2.0;
    double h_2r = container->h - r*2.0;
    double l_2r;
    double bx, by, bz, a;
  
    if(setting->is_3D) {
      l_2r = container->l - r*2.0;

      bx = beta3D_at_x(setting, particleNum, w_2r, h_2r, l_2r, r);
      by = beta3D_at_x(setting, particleNum, h_2r, l_2r, w_2r, r);
      bz = beta3D_at_x(setting, particleNum, l_2r, w_2r, h_2r, r);
      a = alpha3D(particleNum, w_2r, h_2r, l_2r, r);
      //printf("%f %f %f: bx=%f by=%f bz=%f a=%f\n", container->w, container->h, container->l, bx,  by, bz, a);
      analyticalP = (1.0 - (bx + by + bz)/3.0) / (1.0-a);
    } else {
      bx = beta2D_at_x(setting, particleNum, w_2r, h_2r, r);
      by = beta2D_at_x(setting, particleNum, h_2r, w_2r, r);
      a = alpha2D(particleNum, w_2r, h_2r, r);
      //printf("%f %f: bx=%f by=%f a=%f\n", container->w, container->h, bx,  by, a);
      analyticalP = (1.0 - 0.5*(bx + by)) / (1.0-a);
    }
  }

  return analyticalP;
}

void processResult(Setting *setting,
  double *pressureBuf,
  int particleSimulationCount,
  Container *container,
  int particleNum,
  double **inTimeBuf, 
  int positionTimeCount,
  double *ptThicknessBuf,
  double *ptValueBuf,
  double *ptVarBuf
) {
  double pressureAve, pressureStd;
  double analyticalP, analyticalDensity;
  double inTimeAve, inTimeStd;
  int i;

  if(setting->do_density_stat) {
    QuadraticFittingCoef coef;

    for(i=0; i<positionTimeCount; i++) {
      averageAndSd(&inTimeAve, &inTimeStd, 
        inTimeBuf[2*i], particleSimulationCount);
      printf("y2=%f ave=%f var=%f\n", ptThicknessBuf[i], inTimeAve, inTimeStd);
      ptValueBuf[i] = inTimeAve;
      ptVarBuf[i] = inTimeStd*inTimeStd;
    }
    quadraticFitting(&coef, ptThicknessBuf, ptValueBuf,
      ptVarBuf, positionTimeCount);
    printf("density fitting: c0=%f c1=%f c2=%f rmse=%f c0 std dev.=%f\n", 
      coef.c0, coef.c1, coef.c2, coef.rmse, coef.c0Sd);

    analyticalDensity = setting->density_prediction;

    printf("Density: %f Analytical Density: %f Simul-Analytic: %f\n", 
       coef.c0, analyticalDensity, coef.c0-analyticalDensity);

    for(i=0; i<positionTimeCount; i++) {
      averageAndSd(&inTimeAve, &inTimeStd, 
        inTimeBuf[2*i+1], particleSimulationCount);
      printf("y2=%f vNormalSq ave=%f v^2 var=%f\n", 
        ptThicknessBuf[i], inTimeAve, inTimeStd);
      // ptValueBuf[i] = inTimeAve;
    }
    // quadraticFitting(&coef, ptThicknessBuf, ptValueBuf, NULL, 
    //   positionTimeCount);

    // printf("velocity fitting: c0=%f c1=%f c2=%f rmse=%f\n", 
    //   coef.c0, coef.c1, coef.c2, coef.rmse);
  } else {
    analyticalP = analyticalPressure(setting, particleNum, container);
  
    averageAndSd(&pressureAve, &pressureStd, 
      pressureBuf, particleSimulationCount);

    if(setting->is_3D) {
      printf("%f %f %f %f %g %f %g\n", 
       container->w, container->h, container->l, 
       pressureAve, pressureStd, analyticalP, pressureAve - analyticalP);
    } else {
      printf("%f %f %f %g %f %g\n", container->w, container->h, 
       pressureAve, pressureStd, analyticalP, pressureAve - analyticalP);
    }
  }
}

int pressureStatistics(Setting *setting, Container *container, int particleNum) {
  Impulse impulse;
  Particle *p;
  double *pressureBuf;
  double simulationT = SIMULATION_TIME;
  double p0;
  double impulsePerArea;
  int particleSimulationCount = 50;
  int i;

  p = allocateParticles(particleNum);
  pressureBuf = (double *) malloc(sizeof(double)*particleSimulationCount);
  if(!p || !pressureBuf) {
    free(p);
    free(pressureBuf);
    fprintf(stderr, "Failed to allocate particles or pressureBuf\n");
    return 0;
  }

  printf("Particle num: %d, interiorP: %d, "
    "exclude %f each side\n",
    particleNum, setting->is_interior_p,
    setting->excluded_d);
  p0 = pressure_0(setting, container, particleNum);

  for(i=0; i<particleSimulationCount; i++) {
    initImpulse(&impulse);
    oneSimulation(&impulse, setting, p, particleNum, container,
        simulationT, NULL, 0);

    impulsePerArea = impulsePerUnitArea(setting, &impulse, container);
    //printf("impulse %f %f %f\n", impulse.impulseX, impulse.impulseY, impulse.impulseZ);
    //printf("impulse per area %f\n", impulsePerArea);
    pressureBuf[i] = impulsePerArea/simulationT/p0;
  }

  processResult(setting, pressureBuf, particleSimulationCount, container,
      particleNum, NULL, 0, NULL, NULL, NULL);

  free(p);
  free(pressureBuf);

  return 1;
}

int densityStatistics(Setting *setting, Container *container, int particleNum) {
  Particle *p;
  double *pressureBuf;
  double *ptThicknessBuf;
  double *ptValueBuf;
  double *ptVarBuf;
  double **inTimeBuf;
  PositionTime *pt;
  double simulationT = SIMULATION_TIME;
  double p0;
  int particleSimulationCount = 50;
  int positionTimeCount = 8;
  int i, j;

  p = allocateParticles(particleNum);
  pt = allocateDensityStatisticsBoxes(positionTimeCount);
  ptThicknessBuf = (double *) malloc(sizeof(double) * positionTimeCount);
  ptValueBuf = (double *) malloc(sizeof(double) * positionTimeCount);
  ptVarBuf = (double *) malloc(sizeof(double) * positionTimeCount);
  inTimeBuf = allocateInTimeBuf(positionTimeCount, particleSimulationCount);
  if(!pt || !ptThicknessBuf || !ptValueBuf || !ptVarBuf || !inTimeBuf) {
    fprintf(stderr, "Failed to allocate buffers\n");
    free(p);
    free(pt);
    free(ptThicknessBuf);
    free(ptValueBuf);
    free(ptVarBuf);
    freeInTimeBuf(inTimeBuf, 2*positionTimeCount);
    return 0;
  }

  printf("Particle num: %d, "
    "exclude %f each side, "
    "density stat at once: %d\n",
    particleNum, 
    setting->excluded_d,
    setting->density_at_once);
  p0 = pressure_0(setting, container, particleNum);

  initDensityStatisticsBoxes(setting, pt, positionTimeCount, 
    ptThicknessBuf, container, particleNum);

  setting->density_prediction = particleNum > 2 
    ? analyticalBoundaryDensityPrimaryS(setting, particleNum, container)
    : analyticalBoundaryDensityPrimary(setting, particleNum, container);

  printf("Analytical Density Prediction: %f\n", setting->density_prediction);

  if(setting->density_at_once) {
    for(i=0; i<particleSimulationCount; i++) {
      oneSimulation(NULL, setting, p, particleNum, container,
        simulationT, pt, positionTimeCount);

      for(j=0; j<positionTimeCount; j++) {
        inTimeBuf[2*j][i] = pt[j].inTime/
          (simulationT*particleNum*pt[j].volumeRatio);
        inTimeBuf[2*j+1][i] = pt[j].inTimeVySq/pt[j].inTime;
      }
    }
  } else {
    for(j=0; j<positionTimeCount; j++) {
      for(i=0; i<particleSimulationCount; i++) {
        oneSimulation(NULL, setting, p, particleNum, container,
          simulationT, pt+j, 1);

        inTimeBuf[2*j][i] = pt[j].inTime/
          (simulationT*particleNum*pt[j].volumeRatio);
        inTimeBuf[2*j+1][i] = pt[j].inTimeVySq/pt[j].inTime;
      }
    }
  }

  processResult(setting, pressureBuf, particleSimulationCount, container,
    particleNum, inTimeBuf, positionTimeCount, 
    ptThicknessBuf, ptValueBuf, ptVarBuf);

  free(p);
  free(pt);
  free(ptThicknessBuf);
  free(ptValueBuf);
  free(ptVarBuf);
  freeInTimeBuf(inTimeBuf, 2*positionTimeCount); 

  return 1;
}

void print_alpha_beta(Setting *setting, int particleNum, double w, double h, double l, double r) {
  double w_2r = w-2.0*r;
  double h_2r = h-2.0*r;
  double l_2r = l-2.0*r;

  double bx = beta3D_at_x(setting, particleNum, w_2r, h_2r, l_2r, r);
  double by = beta3D_at_x(setting, particleNum, h_2r, l_2r, w_2r, r);
  double bz = beta3D_at_x(setting, particleNum, l_2r, w_2r, h_2r, r);
  double a = alpha3D(particleNum, w_2r, h_2r, l_2r, r);
  printf("%f %f %f %f %f\n", bx, by, bz, a, (1.0-(bx+by+bz)/3.0)/(1.0-a));
}

int main(int argc, char **argv) {
  Container container2D[] = {
    {.w=19, .h=20}, 
    {.w=16, .h=18}, 
    {.w=11, .h=12},
    {.w=19, .h=20}, 
    {.w=16, .h=18}, 
    {.w=19, .h=20}, 
    {.w=16, .h=18}
  };
  Container container3D[] = {
    {.w=18, .h=19, .l=19},
    {.w=14, .h=16, .l=17}, 
    {.w=10, .h=11, .l=12}, 
    {.w=14, .h=16, .l=17}, 
    {.w=10, .h=11, .l=12}, 
    {.w=14, .h=16, .l=17}, 
    {.w=10, .h=11, .l=12}
  };
  int particleNum2D[] = {2, 2, 2, 5, 5, 8, 8};
  int particleNum3D[] = {2, 2, 2, 5, 5, 8, 8};

  Container container2D_density[] = {
    {.w=16, .h=18}, 
    {.w=16, .h=18}, 
    {.w=16, .h=18}
  };
  Container container3D_density[] = {
    {.w=10, .h=11, .l=12}, 
    {.w=10, .h=11, .l=12}, 
    {.w=10, .h=11, .l=12}, 
  };
  int particleNum2D_density[] = {2, 5, 8};
  int particleNum3D_density[] = {2, 5, 8};

  // set random seed so even if we start the run at a middle index,
  // we still get the same result
  int rSeedArr[] = {2134, 4743, 8843, 9932, 7843, 6432, 1473, 6347, 9742, 7438,
     6467, 4242, 4743, 7423, 7538, 7348, 6432, 5372, 6372, 6728};

  Setting setting;
  Container *container;
  int *particleNumArr;
  int particleNumArrCount, particleNum;
  int container2DCount = sizeof(container2D)/sizeof(container2D[0]);
  int container3DCount = sizeof(container3D)/sizeof(container3D[0]);

  initSetting(&setting);
  int pressureStartIndex = 0;
  for(int is_interior_p = 0; is_interior_p < 2; is_interior_p++) {
    setting.is_interior_p = is_interior_p;
    setting.excluded_d = is_interior_p
      ? 2.0 * setting.radius : 0;

    for(int i=0; i<2; i++) {
      if(i==0) {
        setting.is_3D = 0;
        particleNumArr = particleNum2D;
        particleNumArrCount = sizeof(particleNum2D)/sizeof(particleNum2D[0]);
        container = container2D;
      } else {
        setting.is_3D = 1;
        particleNumArr = particleNum3D;
        particleNumArrCount = sizeof(particleNum3D)/sizeof(particleNum3D[0]);
        container = container3D;
      }

      for(int j=pressureStartIndex; j<particleNumArrCount; j++) {
        particleNum = particleNumArr[j];
        if(setting.is_interior_p==0 && particleNum > 2) {
          continue;
        }

        srandom(rSeedArr[j]);
        if(!pressureStatistics(&setting, container+j, particleNum)) {
          return 1;
        }
      }
      printf("\n");
    }
  }

  setting.do_density_stat = 1;
  setting.is_interior_p = 1;
  setting.excluded_d = 2.0*setting.radius;
  int densityStartIdx = 0;
  for(int i=0; i<2; i++) {
    if(i==0) {
      setting.is_3D = 0;
      particleNumArr = particleNum2D_density;
      particleNumArrCount = sizeof(particleNum2D_density)/sizeof(particleNum2D_density[0]);
      container = container2D_density;
    } else {
      setting.is_3D = 1;
      particleNumArr = particleNum3D_density;
      particleNumArrCount = sizeof(particleNum3D_density)/sizeof(particleNum3D_density[0]);
      container = container3D_density;
    }

    for(int j=densityStartIdx; j<particleNumArrCount; j++) {
      particleNum = particleNumArr[j];
      //setting.is_interior_p = particleNum==2 ? 0:1;
      srandom(rSeedArr[j]);
      if(!densityStatistics(&setting, container+j, particleNum)) {
        return 1;
      }
    }
  }

  return 0;
}
