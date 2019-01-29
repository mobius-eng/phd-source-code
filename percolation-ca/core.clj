;; * Core
(ns percolation-ca.core
  (:require [clojure.set :as clj-set]
            [quil.core :as q]
            [quil.middleware :as m]
            [dk.ative.docjure.spreadsheet :as xl]
            [percolation-ca.lattice :as lat]
            [percolation-ca.simulation :as sim]
            [percolation-ca.simulation-visualization :as vis]
            [percolation-ca.plot :as plt]))

;; ** Parameters controlling computation

;; *** Rules
(def ^:dynamic *accept-probability*
  "Probabilities (out of 100) to accept the fluid
  parcel for each type of the pore and for each state"
  {:throat [5 100 100 0]
   :body [20 80 20 0]
   :drainage [20 80 20 0]})

(def ^:dynamic *accept-direction-factor*
  "Direction factor for accepting probabilities. A type of pore
  might be more like to accept the parcel if it comes from a particular direction.
  This is a way to create a preferred direction (for example, down).
  The parameters show the factors FROM a direction"
  {:throat {:right 1.0 :left 1.0 :up 1.0 :down 1.0}
   :body {:right 2.0 :left 2.0 :down 0.1 :up 5.0}
   :drainage {:right 1.0 :left 1.0 :down 0.1 :up 5.0}})

(def ^:dynamic *expell-probability*
  "Probability of expelling a fluid parcel from a particular
  type of pore."
  {:throat   [0 0.02   5  50]
   :body     [0   1  50 100]
   :drainage [0   0   0   0]})

(def ^:dynamic *states*
  "A vector of format
  `index`: number of parcels
  `value`: state number
  If number of parcels exceeds max, the state for the max is assumed"
  [0 1 2 2 2 2 3])

(def ^:dynamic *moving*
  [1 1 1 1 1 1 2 2 5 10 30 100])

;; *** Lattice structure
(def ^:dynamic *m*
  "Number of horizontal layers"
  20)

(def ^:dynamic *n*
  "Number of sites in a horizontal layer (for odd layers: will be -1)"
  25)

(def ^:dynamic *throat-fraction*
  "Fraction of throat pores in the lattice"
  0.7)

;; *** Simulation
(def ^:dynamic *inflow*
  "How much liquid and where goes in"
  (let [n2 1]
    (comment (quot (inc *n*) 2))
    {n2 3, 2 3}))

;; *** Viualization parameters
(def ^:dynamic *frame-rate* 30)

(def ^:dynamic *size*
  "Size of one bond"
  30)

(def ^:dynamic *margin*
  "A margin between visual components"
  10)

;; ** Simulation
;; *** Setup state
(defn setup-state []
  (let [lattice (lat/initialize-lattice *m*
                                        *n*
                                        (fn [_] (< (rand) *throat-fraction*))
                                        :throat)
        rules (sim/make-rules *states*
                              *accept-probability*
                              *accept-direction-factor*
                              *expell-probability*
                              *moving*)
        simulation (vis/setup-state lattice rules *size* *inflow*)]
    (assoc simulation
           :mass-plot {:data [[0 0]]
                       :frame [0 0 300 200]
                       :color [255 0 0]
                       :thickness 2
                       :title "Total number of parcels"
                       :margin *margin*}
           :rtd-plot {:data [[0 0]]
                      :frame [0 0 300 200]
                      :color [0 0 255]
                      :thickness 2
                      :title "RTD of fluid"
                      :margin *margin*}
           :inflow-history [[0 *inflow*]])))

(defn update-state [{:keys [mass-plot rtd-plot] :as model}]
  (let [new-model (vis/update-state model)
        new-mass-plot (update mass-plot
                              :data
                              (fn [data]
                                (conj data
                                      [((new-model :state) :time)
                                       (reduce #(+ %1 (count %2))
                                              0
                                              (vals ((new-model :state)
                                                     :active-bonds)))])))
        new-rtd-plot (assoc rtd-plot
                            :data
                            (cons [0 0]
                             (sort (fn [[x1 _] [x2 _]]
                                     (cond (< x1 x2) -1
                                           (> x1 x2) 1
                                           :else 0))
                                   (seq (-> new-model :state :removed)))))]
    (assoc new-model
           :mass-plot new-mass-plot
           :rtd-plot new-rtd-plot)))

(def state0 (setup-state))
;; (def state1 (reduce (fn [state _] (update-state state)) state0 (range 1000)))

;; state1

;; (seq (-> state0 :rules :accept-probability))

(defn draw [model]
  (q/background 255)
  (vis/draw model)
  (q/with-translation [(+ (* 2 *margin*) (* *n* *size*)) *margin*]
      (plt/draw-plot (model :mass-plot)))
  (q/with-translation [(+ (* 2 *margin*) (* *n* *size*)) (+ (* 2 *margin*) 200)]
    (plt/draw-plot (model :rtd-plot))))

(defn save-state [{:keys [lattice rules state mass-plot rtd-plot inflow-history]}]
  (let [wb (xl/create-workbook "Configuration"
                               [["M" *m*]
                                ["N" *n*]
                                ["Number of bonds" (-> lattice :bonds count)]
                                ["Max bond capacity" (-> rules :states count dec)]
                                ["Throat fraction" *throat-fraction*]]
                               "Rules"
                               (concat
                                [["Accepting"]]
                                (map (fn [[k v]] (cons (str k) v))
                                     (rules :accept-probability))
                                [["Expelling"]]
                                (map (fn [[k v]] (cons (str k) v))
                                     (rules :expell-probability))
                                [["Moving degeneration factor"]]
                                [(rules :moving)])
                               "Inflow History"
                               (cons ["Generation" "Inflow"]
                                     (map (fn [[t inflow]]
                                            [t (reduce (fn [v [_ r]]
                                                         (+ v r))
                                                       0
                                                       inflow)])
                                          inflow-history))
                               "Mass"
                               (cons ["Generation" "Parcels"]
                                     (mass-plot :data))
                               "RTD"
                               (cons ["Generations spent" "Parcels"]
                                     (rtd-plot :data)))]
    (xl/save-workbook!
     (str "results_"
          (.format (new java.text.SimpleDateFormat "yyyy-MM-dd_HH-mm-ss") (java.util.Date.))
          ".xlsx")
     wb)))

(defn handle-key [state event]
  (case (event :raw-key)
    \s (-> state
           (update :inflow-history #(conj % [(-> state :state :time) []]))
           (assoc :inflow []))
    \f (-> state
           (update :inflow-history #(conj % [(-> state :state :time) *inflow*]))
           (assoc :inflow *inflow*))
    \r (assoc-in state [:state :removed] {})
    \R (assoc-in state [:mass-plot :data] [])
    \C (-> state (assoc-in [:state :active-bonds] {}) (assoc-in [:state :parcels] #{}))
    \K (setup-state)
    \p (do (q/save-frame (str "percolation-"
                              (.format (new java.text.SimpleDateFormat
                                            "yyyy-MM-dd_HH-mm-ss")
                                       (java.util.Date.))
                              ".png"))
           state)
    \x (do (save-state state) state)
    \u (assoc state :rules (sim/make-rules *states*
                                           *accept-probability*
                                           *accept-direction-factor*
                                           *expell-probability*))
    state))

(defn main []
  (q/defsketch percolation-ca
    :size [(+ (+ 20 (* *n* *size*)) (* 3 *margin*) 300)
           (+ 20 (* *m* (* *size* (Math/sin (/ Math/PI 3)))))]
    :setup (fn[]
             (q/frame-rate *frame-rate*)
             ;; (q/text-font (q/create-font "Arial" 12 true))
             (setup-state))
    :update (comp update-state update-state)
    :draw draw
    :middleware [m/fun-mode]
    :key-typed handle-key))

(main)
