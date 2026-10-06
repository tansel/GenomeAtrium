/*
 * Remote calls that survive a flaky connection. retry(fn) runs fn (a function returning a
 * promise) up to `tries` times, waiting 1, 2, 4... seconds between attempts, and rejects with
 * the last error. Used for the libraries (three.js, GenomeSpy) and the protein lookups, where
 * one dropped request used to leave a view broken until the page was reloaded.
 */
(function (G) {
  function retry(fn, tries, firstDelay) {
    tries = tries || 4; var delay = firstDelay || 1000;
    return new Promise(function (resolve, reject) {
      (function attempt(n) {
        Promise.resolve().then(fn).then(resolve, function (err) {
          if (n + 1 >= tries) { reject(err); return; }
          setTimeout(function () { attempt(n + 1); }, delay * Math.pow(2, n));
        });
      })(0);
    });
  }
  // fetch that treats an HTTP 429 or 5xx as a failure worth retrying (4xx other than 429 is final)
  function fetchRetry(url, opts, tries) {
    return retry(function () {
      return fetch(url, opts).then(function (r) {
        if (r.status === 429 || r.status >= 500) throw new Error('HTTP ' + r.status + ' from ' + url.split('/')[2]);
        return r;
      });
    }, tries);
  }
  G.net = { retry: retry, fetchRetry: fetchRetry };
})(globalThis.G = globalThis.G || {});
