(ns percolation-ca.plot
  (:require [quil.core :as q]
            [quil.middleware :as m]))

(def ^:dynamic *title-space* 20)

(defn- screen-coordinates
  "Return screen coordinates of the data in the plot"
  [{:keys [data frame margin]}]
  (let [[left top width height] frame
        [xmin ymin] (if (empty? data) [0 0] (first data))
        [xmin xmax ymin ymax] (reduce (fn [[cxmin cxmax cymin cymax] [x y]]
                                        [(min cxmin x)
                                         (max cxmax x)
                                         (min cymin y)
                                         (max cymax y)])
                                      [xmin xmin ymin ymin]
                                      data)
        xrange (max (- xmax xmin) 0.0001)
        yrange (max (- ymax ymin) 0.0001)
        true-height (- height *title-space*)]
    (map (fn [[x y]]
           [(+ margin left (q/round (* (/ (- width margin margin) xrange) (- x xmin))))
            (+ (- margin) top height (q/round (* (/ (- true-height margin margin) yrange) (- ymin y))))])
         data)))

(defn draw-plot [{:keys [data frame title color thickness margin] :as plot}]
  (let [[left top width height] frame
        screen-data (partition 2 1 (screen-coordinates plot))]
    (q/stroke-weight 1)
    (q/with-stroke [5 5 5]
      (q/rect left top width height))
    (q/with-fill [0 0 0]
      (q/text title (+ left margin) (+ top *title-space*)))
    (q/stroke-weight thickness)
    (q/with-stroke color
      (doseq [[[x1 y1] [x2 y2]] screen-data]
        (q/line x1 y1 x2 y2)))))


(defn setup-example []
  {:data [[0.0 0.0]]
   :frame [0 0 300 200]
   :title "Sin(x) example"
   :color [255 0 0]
   :thickness 3
   :margin 5})

(defn update-example [{:keys [data] :as plot}]
  (let [[x y] (last data)]
    (update plot :data (fn [data]
                         (conj data [(+ x 0.05)
                                     (Math/sin (+ x 0.05))])))))

(defn run-example []
  (q/defsketch plot-example
    :size [600 600]
    :setup setup-example
    :update update-example
    :draw draw-plot
    :middleware [m/fun-mode]))

;; (run-example)
