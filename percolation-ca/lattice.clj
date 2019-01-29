(ns percolation-ca.lattice
  (:require [clojure.set :as clj-set]))

;; * Lattice structure
;; ** Bonds
(defn bond-neighbours
  "Returns all the bonds that have common sites with the specified bond"
  [bond-index lattice]
  (((lattice :bonds) bond-index) :neighbours))


(defn- specialize-bond-map
  "INTERNAL
  Bond map is the map of the form
  {:left l :right r :up-left ul :up-right ur :down-left dl :down-right dr}
  This function trims the map of non-existent keys for particular sites"
  [bond-map m n i j]
  (let [jmax-even (dec n)
        jmax-odd (dec jmax-even)
        imax (dec m)
        bond-map! (transient bond-map)]
    (persistent!
     (if (even? i)
       (let [bond-map! (cond (= j 0) (dissoc! bond-map!
                                              :left :up-left :down-left)
                             (= j jmax-even) (dissoc! bond-map!
                                                      :right :up-right :down-right)
                             :else bond-map!)]
         (cond (= i 0) (dissoc! bond-map! :up-left :up-right)
               (= i imax) (dissoc! bond-map! :down-left :down-right)
               :else bond-map!))
       (let [bond-map! (cond (= j 0) (dissoc! bond-map! :left)
                             (= j jmax-odd) (dissoc! bond-map! :right)
                             :else bond-map!)]
         (cond (= i imax) (dissoc! bond-map! :down-left :down-right)
               :else bond-map!))))))

(defn- bonds-add-site
  "INTERNAL
  Add site number `siten` as a site to bonds specified by numbers `nbonds`"
  [bonds siten bond-map]
  (reduce (fn [acc val]
            (update acc val 
                    (fn [bond] (update bond :sites #(conj % siten)))))
          bonds
          (vals bond-map)))

(defn- site-add-bonds
  "INTERNAL
  Add bonds specified by the bond-map to the site numbered `siten`"
  [sites siten bonds-map]
  (update sites siten (fn [site] (update site :bonds #(merge % bonds-map)))))

(defn- make-bond-map
  "INTERNAL
  Creates a general bond-map for the site. WARNING: this map will containt
  non-existing bonds (for boundary sites), these bonds need to be trimmed
  using `specialize-bond-map`"
  [n llb i j]
  (let [nbond-left (+ llb j)
        nbond-right (inc nbond-left)
        nbond-up-left (+ (- llb (* 2 (dec n))) (* 2 j) (if (odd? i) 1 0))
        nbond-up-right (inc nbond-up-left)
        nbond-down-left (+ llb n -1 (* 2 j))
        nbond-down-right (inc nbond-down-left)]
    {:left nbond-left
     :right nbond-right
     :up-left nbond-up-left
     :up-right nbond-up-right
     :down-left nbond-down-left
     :down-right nbond-down-right}))

(def bond-directions
  "INTERNAL
  Represents the direction of bonds in vector form (used for convenience)"
  {:up-left [-1 -1] :up-right [-1 +1]
   :left [0 -1] :right [0 +1]
   :down-left [+1 -1] :down-right [+1 +1]})

(defn make-bonds
  "LOW-LEVEL
  Creates trivial bonds for the rectangular lattice.
  All the bonds are by default of type `:body`. No information on sites
  is added yet. For a higher level generation of the lattice see e.g.
  `make-lattice`"
  [m n]
  (let [n-bonds (+ (* (- (* 3 n) 3) (quot m 2))
                   (*(- (* 3 n) 4) (quot (dec m) 2))
                   n -2
                   (if (even? (dec m)) 1 0))]
    (mapv (fn [n]
            {:index n
             :sites #{}
             :neighbours {}
             :type :body})
          (range n-bonds))))

;; ** Sites
(defn make-sites
  "LOW-LEVEL
  Generates the sites for the rectangular lattice. Sites are not fully initialized:
  bonds are not set. Site's index and lattice location are set."
  [m n]
  (let [n-sites (- (* m n) (quot m 2))
        imax (dec m)
        jmax-even (dec n)
        jmax-odd (dec jmax-even)]
    (loop [n 0 i 0 j 0 sites (transient [])]
      (if (= n n-sites)
        (persistent! sites)
        (if (or (and (even? i) (= j jmax-even)) (and (odd? i) (= j jmax-odd)))
          (recur (inc n) (inc i) 0 (conj! sites
                                          {:index n
                                           :loc [i j]
                                           :bonds {}}))
          (recur (inc n) i (inc j) (conj! sites
                                          {:index n
                                           :loc [i j]
                                           :bonds {}})))))))


(defn- bond-neighbour-direction [[v h]]
  (cond (neg? v) :up
        (pos? v) :down
        (neg? h) :left
        :else :right))

(defn- negate [[x y]] [(- x) (- y)])
(defn- add [[x y] [p q]] [(+ x p) (+ y q)])


;; TODO: busy fixing
(defn- add-bond-neighbours
  "Adds bond neighbours to the bond. Neighbours is a map
  of bond index to one of the directions :up, :down, :left, :right"
  [bond-index bonds lattice-sites]
  (let [bond (bonds bond-index)
        sites (bond :sites)
        neighbours (reduce (fn [m s]
                             (let [site-bonds ((lattice-sites s) :bonds)
                                   site-dir (->> site-bonds
                                                 (filter #(= (val %) bond-index))
                                                 first
                                                 key
                                                 bond-directions
                                                 negate)]
                               (reduce (fn [m [bond-dir site-bond]]
                                         (if (not= bond-index site-bond)
                                           (let [neighbour-dir (->> bond-dir
                                                                    bond-directions
                                                                    (add site-dir)
                                                                    bond-neighbour-direction)]
                                             (assoc m site-bond neighbour-dir))
                                           m))
                                       m
                                       site-bonds)))
                           {}
                           sites)]
    (assoc-in bonds [bond-index :neighbours] neighbours)))

;; ** Structured lattice

;; Lattice itself is a vector of bonds. Each bond is numbered. These
;; numbers need to be translated to bond position, i.e. location of bond
;; sites.
(defn- make-lattice-layer [m n i siten llb bonds sites]
  (let [llb-inc (- (* 3 (dec n)) (if (odd? i) 1 0))
        jmax (if (even? i) (dec n) (- n 2))]
    (loop [j 0 siten siten bonds bonds sites sites]
      (let [bond-map (-> (make-bond-map n llb i j)
                         (specialize-bond-map m n i j))]
        (cond (> j jmax) [siten (+ llb llb-inc) bonds sites]
              :else (recur (inc j)
                           (inc siten)
                           (bonds-add-site bonds siten bond-map)
                           (site-add-bonds sites siten bond-map)))))))

(defn make-lattice
  "Constructs fully initialized rectangular lattice with `m` sites
  in vertical direction and `n` sites in horisontal direction (due to
  triangular nature of the lattice, the number of horisontal sites
  is `n` for all even horisontal layers[e.g., 0, 2, ...] and `n-1` for
  all odd layers [e.g. 1,3,...])"
  [m n]
  (loop [i 0 siten 0 llb -1 bonds (make-bonds m n) sites (make-sites m n)]
    (cond (= i m) {:size [m n]
                   :bonds (reduce (fn [new-bonds bond]
                                    (add-bond-neighbours (bond :index) new-bonds sites))
                                  bonds
                                  bonds)
                   :sites sites}
          :else (let [[siten llb bonds sites] (make-lattice-layer m n
                                                                  i
                                                                  siten
                                                                  llb
                                                                  bonds
                                                                  sites)]
                  (recur (inc i) siten llb bonds sites)))))

(defn bottom-bonds
  "Returns the range of bottom bonds indices"
  [lattice]
  (let [[m n] (lattice :size)
        bmax (count (lattice :bonds))
        bottom (if (even? m) (- n 2) (- n 1))]
    (range (- bmax bottom) bmax)))

;; ** Generate throats
(defn update-bonds-type [lattice update? new-type]
  (update lattice
          :bonds
          (fn [bonds]
            (mapv (fn [bond]
                    (if (update? bond)
                      (assoc bond :type new-type)
                      bond))
                  bonds))))

(defn make-drainage-bonds [{:keys [size bonds] :as lattice}]
  (let [bmax (count bonds)
        [m n] size
        bottom-bonds (if (even? m) (- n 2) (- n 1))]
    (update lattice
            :bonds
            (fn [bonds]
              (loop [indices (range (- bmax bottom-bonds) bmax)
                     bonds (transient bonds)]
                (if (empty? indices)
                  (persistent! bonds)
                  (recur
                   (rest indices)
                   (let [index (first indices)]
                     (assoc! bonds
                             index
                             (assoc (bonds index) :type :drainage))))))))))

(defn initialize-lattice [m n update-bonds-type? new-type]
  (-> (make-lattice m n)
      (update-bonds-type update-bonds-type? new-type)
      make-drainage-bonds))

