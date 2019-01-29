(ns percolation-ca.simulation-visualization
  (:require [percolation-ca.lattice :as lat]
            [percolation-ca.simulation :as sim]
            [quil.core :as q]
            [quil.middleware :as m]))

;; *** Uitilities

(defn calculate-site-points [sites bond-size]
  (let [layer-height (* bond-size (Math/sin (/ Math/PI 3)))
        odd-i-shift (* bond-size (Math/cos (/ Math/PI 3)))]
   (mapv (fn [site]
           (let [[i j] (site :loc)]
             [(+ 10 (* j bond-size) (if (odd? i) odd-i-shift 0))
              (+ 10 (* i layer-height))]))
         sites)))

;; (max 0 (- 200 (* parcels 40)))
(defn draw-bond [{:keys [lattice state site-points]} bond]
  (let [[s1 s2] (seq (((lattice :bonds) bond) :sites))
        [x1 y1] (site-points s1)
        [x2 y2] (site-points s2)
        bond-type (((lattice :bonds) bond) :type)
        parcels (count ((state :active-bonds) bond))
        color (max 0 (- 180 (* parcels 30)))]
    (q/with-stroke (case bond-type
                     :body [color 255 color]
                     :throat [color color 255]
                     :drainage [255 color color]
                     [50 50 50])
      (q/line x1 y1 x2 y2))))

(defn draw [{:keys [lattice state removed] :as model}]
  (q/stroke-weight 1)
  (doseq [bond (state :active-bonds)]
    (draw-bond model (key bond)))
  (let [[body-parcels throat-parcels] (reduce (fn [[b t] [bond parcels]]
                                                (case (((lattice :bonds) bond) :type)
                                                  :body [(+ b (count parcels)) t]
                                                  :throat [b (+ t (count parcels))]
                                                  [b t]))
                                              [0 0]
                                              (state :active-bonds))]
    (q/with-stroke [0 0 0 0]
      (q/with-fill [150 150 150 100]
        (q/rect 5 5 150 160)))
    (q/with-fill [0 0 0]
        (q/text (str "Generation: " (state :time)) 10 30)
        (q/text (str "Body parcels: " body-parcels) 10 60)
        (q/text (str "Throat parcels: " throat-parcels) 10 90)
        (q/text (str "Frame rate: " (q/round (q/current-frame-rate))) 10 150))))

(defn setup-state [lattice rules bond-size inflow]
  (let [model (sim/make-simulation lattice rules)
        site-points (calculate-site-points (lattice :sites) bond-size)]
    (assoc model
           :site-points site-points
           :inflow inflow)))

(defn update-state [{:keys [lattice state inflow] :as model}]
  (let [new-model (sim/time-step model inflow)]
   new-model))

