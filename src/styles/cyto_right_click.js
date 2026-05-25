(function () {
    function isInsideCytoscape(target) {
        const cyto = document.getElementById("reactome-network");
        return cyto && cyto.contains(target);
    }

    document.addEventListener(
        "contextmenu",
        function (event) {
            if (!isInsideCytoscape(event.target)) {
                return;
            }

            event.preventDefault();

            if (
                !window.dash_clientside ||
                !window.dash_clientside.set_props
            ) {
                console.warn("dash_clientside.set_props is not available.");
                return;
            }

            console.log("RIGHT CLICK CAPTURED", event.clientX, event.clientY);

            window.dash_clientside.set_props("right-click-event-store", {
                data: {
                    x: event.clientX,
                    y: event.clientY,
                    ts: Date.now(),
                },
            });
        },
        true
    );
})();