#lang racket
(require "declarations.rkt")
(provide buildTree calcForces moveparticles)

(define (buildTree initialArea particles)
  (let* [(xi (bbox-llx initialArea))
         (yi (bbox-lly initialArea))
         (xf (bbox-rux initialArea))
         (yf (bbox-ruy initialArea))
         (mid-x (/ (+ xi xf) 2))
         (mid-y (/ (+ yi yf) 2))
         (box1 (bbox xi mid-y mid-x yf))
         (box2 (bbox mid-x mid-y xf yf))
         (box3 (bbox xi yi mid-x mid-y))
         (box4 (bbox mid-x yi xf mid-y))
         (particles-in-box1 (filter (lambda (a) (and (< (x-of a) mid-x) (>= (y-of a) mid-y))) particles))
         (particles-in-box2 (filter (lambda (a) (and (>= (x-of a) mid-x) (>= (y-of a) mid-y))) particles))
         (particles-in-box3 (filter (lambda (a) (and (< (x-of a) mid-x) (< (y-of a) mid-y))) particles))
         (particles-in-box4 (filter (lambda (a) (and (>= (x-of a) mid-x) (< (y-of a) mid-y))) particles))
         (mass-in-whole-box (mass-in-box particles))
         (centroid-of-whole-box (centroid-of-particles particles))]
    (cond [(null? particles) (gnode 0 (vec 0 0) '())]
          [(singleton particles) (gnode (particle-mass (car particles)) (particle-posn (car particles)) '())]
          [else (gnode mass-in-whole-box centroid-of-whole-box
                       (list (buildTree box1 particles-in-box1)
                             (buildTree box2 particles-in-box2)
                             (buildTree box3 particles-in-box3)
                             (buildTree box4 particles-in-box4)))])))

(define (x-of part)
  (let* [(position (particle-posn part))]
    (vec-x position)))
(define (y-of part)
  (let* [(position (particle-posn part))]
    (vec-y position)))

(define (mass-in-box particles)
  (cond [(null? particles) 0]
        [else (+ (particle-mass (car particles)) (mass-in-box (cdr particles)))]))

(define (centroid-of-particles particles)
  (if (null? particles) (vec 0 0)
  (let* [(x-sum (sum-of-x particles))
         (y-sum (sum-of-y particles))
         (total-mass (mass-in-box particles))]
    (vec (/ x-sum total-mass) (/ y-sum total-mass)))))
(define (sum-of-x particles)
  (let [(x-sums (map (lambda (a) (* (x-of a) (particle-mass a))) particles))]
    (foldr + 0 x-sums)))
(define (sum-of-y particles)
  (let [(y-sums (map (lambda (a) (* (y-of a) (particle-mass a))) particles))]
    (foldr + 0 y-sums)))

(define (calcForces initialArea tree particles)
  (let* [(l (- (bbox-rux initialArea) (bbox-llx initialArea)))]
    (map (lambda (a) (forces tree l a)) particles)))
(define (forces tree l part)
  (let* [(posn-of-particle (particle-posn part))
         (posn-of-tree (gnode-posn tree))
         (d (distance-between posn-of-tree posn-of-particle))
         (m1 (particle-mass part))
         (m2 (gnode-mass tree))
         (x-component (- (vec-x posn-of-tree) (vec-x posn-of-particle)))
         (y-component (- (vec-y posn-of-tree) (vec-y posn-of-particle)))
         (net-force (if (= d 0) 0 (/ (* g m1 m2) (* d d d))))
         (subtrees (gnode-subtrees tree))]
    (cond [(= 0 d) (vec 0 0)]
          [(> (/ d l) theta) (vec (* net-force x-component) (* net-force y-component))]
          [(null? subtrees) (vec (* net-force x-component) (* net-force y-component))]
          [else (add-forces (forces (car subtrees) (/ l 2) part)
                            (forces (cadr subtrees) (/ l 2) part)
                            (forces (caddr subtrees) (/ l 2) part)
                            (forces (cadddr subtrees) (/ l 2) part))])))

(define (add-forces vec1 vec2 vec3 vec4)
  (let* [(new-vec-x (+ (vec-x vec1) (vec-x vec2) (vec-x vec3) (vec-x vec4)))
         (new-vec-y (+ (vec-y vec1) (vec-y vec2) (vec-y vec3) (vec-y vec4)))]
    (vec new-vec-x new-vec-y)))

(define (distance-between point1 point2)
  (let* [(x1 (vec-x point1))
         (x2 (vec-x point2))
         (y1 (vec-y point1))
         (y2 (vec-y point2))]
    (sqrt (+ (* (- x1 x2) (- x1 x2)) (* (- y1 y2) (- y1 y2))))))

(define (moveparticles particles forces)
  (zipwith new-posn particles forces))
(define (new-posn part force)
  (let* [(x (x-of part))
         (y (y-of part))
         (force-x (vec-x force))
         (force-y (vec-y force))
         (vel-x (vec-x (particle-velocity part)))
         (vel-y (vec-y (particle-velocity part)))
         (mass (particle-mass part))
         (acc-x (/ force-x mass))
         (acc-y (/ force-y mass))
         (new-x (+ x (* vel-x timeslice) (* (/ 1 2) acc-x (* timeslice timeslice))))
         (new-y (+ y (* vel-y timeslice) (* (/ 1 2) acc-y (* timeslice timeslice))))
         (new-vel-x (+ (* acc-x timeslice) vel-x))
         (new-vel-y (+ (* acc-y timeslice) vel-y))]
    (particle mass (vec new-x new-y) (vec new-vel-x new-vel-y))))
