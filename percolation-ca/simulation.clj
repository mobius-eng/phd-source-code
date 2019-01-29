(ns percolation-ca.simulation
  (:require [percolation-ca.lattice :as lat]
            [clojure.set :as clj-set]
            ;; [proto-repl.saved-values :as saved]
            ))

;; * Percolation Cellular Automate

;; ** Rules
(defn make-rules
  "Constructs rules.
  `states` is the vector mapping number of parcels to the state number.
  `accept-probability` is the map of pore types to the vector of accepting probabilities
  `accept-direction-factor` is the map of pore types to the map of direction factors adjusting probabilities
  `expell-probability` is the map of pore types to the vector of expelling probabilities
  `moving` is the vector of probability adjustments for moving parcel"
  [states accept-probability accept-direction-factor expell-probability moving]
  {:states states
   :accept-probability accept-probability
   :accept-direction-factor accept-direction-factor
   :expell-probability expell-probability
   :moving moving})

;; ** Parcels
(defn make-parcel
  "Creates a fluid parcel entering bond `bond-index`"
  [time-in ibond]
  {:id (rand-int 1000000)
   :time-in time-in
   :bond-index ibond})

(defn move-parcel [parcel new-bond-index]
  (assoc parcel :bond-index new-bond-index))

;; ** State of simulation
(defn make-empty-state []
  {:time 0 :parcels #{} :active-bonds {} :removed {}})

(defn make-simulation [lattice rules]
  {:lattice lattice
   :rules rules
   :state (make-empty-state)})

;; ** Get the bond state
(defn bond-state
  "Returns the state of the bond based on parcels present"
  [rules parcels]
  (get-in rules [:states (count parcels)] (last (rules :states))))

;; ** Calculate movement probabilities
(defn accept-probability
  "Calculates the probability of the bond `to` accepting the fluid parcel
  from bond `from` on the lattice `lattice`"
  [state lattice rules ito ifrom]
  (let [to ((lattice :bonds) ito)
        from ((lattice :bonds) ifrom)
        p (get-in rules
                  [:accept-probability
                   (to :type)
                   (bond-state rules (get (state :active-bonds) ito #{}))]
                  0)
        dir-factor (get-in rules [:accept-direction-factor
                                  (to :type)
                                  (get-in to [:neighbours ifrom])])]
    (* p dir-factor)))

(defn expell-probability
  "Calculates the probability of the bond to expell fluid parcel"
  [state lattice rules ifrom]
  (let [bs (bond-state rules (get-in state [:active-bonds ifrom]))])
  (get-in rules
          [:expell-probability
           (get-in lattice [:bonds ifrom :type])
           (bond-state rules (get-in state [:active-bonds ifrom]))]
          100))

(defn expelling?
  "Based on the bond state, probabilitistically decides
  if the bond will be expelling the fluid parcel"
  [expell-probability]
  (< (rand) (/ expell-probability 100)))

;; ** Time update
(declare
 wrandi moved-probability-reduction
 time-step-parcel-update
 apply-inflow apply-drainage)

(defn time-step
 "One time step of simulation. Produces a new simulation
  with the updated state.
  Accepts `simulation` and the map `inflow`, providing the correspondance
  of bond-index to number of parcels entering at the particular time step"
 [{:keys [lattice state rules] :as simulation} inflow]
 (let [prepared-state (-> state
                          (apply-drainage lattice)
                          (apply-inflow lattice inflow))
       new-state! (reduce #(time-step-parcel-update %1 prepared-state lattice rules %2)
                           (-> prepared-state
                              (update :time inc)
                              (assoc :active-bonds (transient {}))
                              (assoc :parcels (transient #{})))
                          (prepared-state :parcels))]
  (assoc simulation :state (-> new-state!
                               (update :active-bonds persistent!)
                               (update :parcels persistent!)))))

(defn wrandi
 "**LOW-LEVEL UTILITY**
  Weighted random index. Returns nil is all weights are zeros"
 [weights]
 (let [total (apply + weights)
       r (* total (rand))]
   (loop [i 0 sum 0 w weights]
     (if (empty? w)
       nil
       (let [sum (+ sum (first w))]
         (if (< r sum)
           i
           (recur (inc i) sum (next w))))))))

(defn- moved-probability-reduction [{:keys [moving]} p nmoved]
  (/ p (get moving nmoved 10000)))

(defn- time-step-parcel-update
  "LOW-LEVEL
   Tracks the changes of the parcel over one time step.
   Incorporates the changes into `updates`"
  [updates state lattice rules {:keys [bond-index] :as parcel}]
  (letfn [(incorporate [u p]
            (-> u
                (update :parcels #(conj! % p))
                (update :active-bonds #(assoc! %
                                               (p :bond-index)
                                               (conj (get-in u [:active-bonds (p :bond-index)] #{})
                                                     p)))))]

   (loop [bond ((lattice :bonds) bond-index) nmoved 0 parcel parcel]
     (if (expelling? (moved-probability-reduction
                      rules
                      (expell-probability state lattice rules (bond :index))
                      nmoved))
       (let [neighbours (keys (lat/bond-neighbours (bond :index) lattice))
             probabilities (map #(accept-probability state lattice rules % (bond :index)) neighbours)]
         (if-let [ind (wrandi probabilities)]
           (let [accepting (nth neighbours ind)]
             (recur ((lattice :bonds) accepting) (inc nmoved) (move-parcel parcel accepting)))
           (incorporate updates parcel)))
       (do
         (when (> nmoved 4) (println "Moved " nmoved " times"))
         (incorporate updates parcel))))))

(defn- introduce-parcel [[parcels! bonds!] time-in index _]
  (let [p (make-parcel time-in index)]
    [(conj! parcels! p)
     (assoc! bonds! index (conj (or (bonds! index #{})) p))]))

;; FIXME
(defn apply-inflow [state lattice inflow]
 (let [time-in (state :time)
       [new-parcels new-bonds] (reduce (fn [pb! [index amount]]
                                        (reduce #(introduce-parcel %1 time-in index %2)
                                                pb!
                                                (range amount)))
                                       [(transient #{}) (transient {})]
                                       inflow)]
  (-> state
      (update :active-bonds #(merge-with clj-set/union % (persistent! new-bonds)))
      (update :parcels clj-set/union (persistent! new-parcels)))))

(defn apply-drainage [state lattice]
  (let [bottom-bonds (lat/bottom-bonds lattice)
        current-time (state :time)]
    (reduce (fn [{:keys [active-bonds] :as state} index]
              (if-let [bond-parcels (active-bonds index)]
                (-> state
                    (update :active-bonds #(dissoc % index))
                    (update :parcels #(clj-set/difference % bond-parcels))
                    (update :removed
                            (fn [removed]
                             (reduce (fn [rems {:keys [time-in]}]
                                      (let [time-spent (- current-time time-in)]
                                       (update rems time-spent #(inc (or % 0)))))
                                     removed
                                     bond-parcels))))
                state))
            state
            bottom-bonds)))
